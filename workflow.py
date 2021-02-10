#! /usr/bin/env python3
from gwf import *
import glob
import os 
#from cutils import *
import json

import pandas as pd
import yaml
import subprocess, sys



gwf = Workflow(defaults={
    #"mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
    "account": "clinicalmicrobio",
    "memory": '1g',
    "walltime": "01:00:00"
})


print("""
                        _            _  ___  
                  _ __ (_)_ __   ___/ |/ _ \\ 
                 | '_ \\| | '_ \\ / _ \\ | (_) |
                 | |_) | | |_) |  __/ |\\__, |
                 | .__/|_| .__/ \\___|_|  /_/ 
~~~~~~~~~~~~~~~~~|_|~~~~~|_|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

""")

"""



"""

# Read config
with open('config.yml') as file:
    config = yaml.full_load(file)

input_base = config['input_base']
output_base = config['output_base']



# Read list of input paths
df = pd.read_table(config['input_list_file'], sep="\t", names = ["batch", "path", "format_specifier"], comment = "#")

# Read list of batches which are completed and hence can be ignored
batches_done = pd.read_table(config['batches_done_file'], sep = "\t", names = ["batch"])['batch'].tolist()


# Prints a nice ending for each joblog.
default_end = f"""echo; echo JOBID $SLURM_JOBID; jobinfo $SLURM_JOBID; echo OK"""


# This is a routine function for 
def conda(env):
    return f"source {config['conda_base']}/etc/profile.d/conda.sh; conda activate {config['conda_env']}"



print("input list file:")
print(df)
print("//")
print()

# TODO: Solve the problem, that this job could theoretically run simultaneously with the nextclade/pangolin targets.



# Iterate over each line in config['input_list_file']
for index, row in df.iterrows():


    print(f"Batch {index}: {row['batch']}", end = "")

    if row["batch"] in batches_done:
        print(f" marked done, ignoring.")
        continue

    print(' not done. Generating jobs:')
    print()
    print(f" path: {row['path']}")


    # Når disse format_specifiers opdateres, er det vigtigt at 
    if row["format_specifier"] == "formatA" or row["format_specifier"] == "formatB": 
        mads_year = 21
    elif row["format_specifier"] == "formatC" or row["format_specifier"] == "formatD": 
        mads_year = 20
    elif row["format_specifier"] == "formatE":
        print(f"\nwarning: refusing to run the batch input list r-script (scripts/parse_path.r), since the input is from mixed years.")
        mads_year = 70 # 70 because it is the least probable year to observe a sample from (backwards and forwards in time)
    else:
        raise Exception(f"format specifier: {row['format_specifier']} is not supported.")

    batch_long = f"{row['batch']}" # e.g. 210108

    print()

    batch_input_file = f"{input_base}/{row['batch']}.tab"

    # Check if the batch_input_file has already been written
    if os.path.exists(batch_input_file):
        print(f" {batch_input_file} has already been written")
    else:
        print(f" creating the file {batch_input_file}")
        try:
            command = f"singularity run {config['singularity_images']}/tidyverse_latest.sif Rscript scripts/parse_path.r {row['batch']} {row['path']} {mads_year} FALSE {row['format_specifier']} > other/input_tmp.tab && mv other/input_tmp.tab {batch_input_file}" # TODO: delete the batch_input_file if it is empty 
            subprocess.run(command, shell = True, check = True)
        except subprocess.CalledProcessError as e:
            print(f"\nAn error occured while initializing {row['batch']}:\n", e)
            sys.exit()


    # Now the batch_input_file file has surely been written, and we can start the actual pipeline
    batch_df = pd.read_table(batch_input_file, sep = "\t", dtype = str)

    print()
    print(batch_df)
    print("//")
    print()



    # TODO: Check that the batch given in input list is similar to the batch specified in the batch_input_file
    batch_sample_list = []
    for batch_index, batch_row in batch_df.iterrows():

        sample_name = f"{batch_row['raw_sample_name']}"
        full_name = f"{row['batch']}.{batch_row['plate']}.{batch_row['moma_serial']}_{batch_row['raw_sample_name']}"
        full_name_clean = full_name.replace(".", "_") # Because gwf or slurm somehow is not compatible with dots!?

        print(full_name, end = " ")



        #print(" ", [batch_row['path'] + i for i in (batch_row['R1'] + batch_row['R2']).split(" ")]); exit()
        t_cat = gwf.target(f"cat__{full_name_clean}",
            inputs = {'forward': [batch_row['path'] + i for i in batch_row['R1'].split(" ")],
                      'reverse': [batch_row['path'] + i for i in batch_row['R2'].split(" ")]}, 
            outputs = {'dir': f"{output_base}/{full_name}/trimmed_reads/",
                       'files': [f"{output_base}/{full_name}/trimmed_reads/{full_name}_val_1.fq.gz",
                                 f"{output_base}/{full_name}/trimmed_reads/{full_name}_val_2.fq.gz"]},
            cores = 4)
        t_cat << \
            f"""

            {conda(config['conda_env'])}
            

            mkdir -p {t_cat.outputs['dir']}
            

            tmp_forward="{output_base}/{full_name}/trimmed_reads/{full_name}_R1{batch_row['extension']}"
            tmp_reverse="{output_base}/{full_name}/trimmed_reads/{full_name}_R2{batch_row['extension']}"


            # Cat the reads together
            cat {" ".join(t_cat.inputs['forward'])} > $tmp_forward
            cat {" ".join(t_cat.inputs['reverse'])} > $tmp_reverse


            # Trim the reads
            trim_galore --paired --fastqc --cores 4 --gzip -o {t_cat.outputs['dir']} --basename {full_name} $tmp_forward $tmp_reverse


            # TODO: Consider removing the catted reads (tmp_).
            # rm $tmp_forward
            # rm $tmp_reverse


            {default_end}

            """

       
        # Map 
        t_map = gwf.target(f"map__{full_name_clean}",
            inputs = t_cat.outputs['files'],
            outputs = {'dir': f"{output_base}/{full_name}/aligned",
                       'bam': f"{output_base}/{full_name}/aligned/{full_name}.sorted.trimmed.bam"},
            cores = 4,
            memory = '4gb',
            walltime = '03:00:00')
        t_map << \
            f"""

            {conda(config['conda_env'])} # TODO: remove bwa from pipe19_a

            mkdir -p {t_map.outputs['dir']}


            mapped="{output_base}/{full_name}/aligned/{full_name}.sorted.tmp.bam"
            renamed="{output_base}/{full_name}/aligned/{full_name}.sorted.bam"
            tmptrimmed="{output_base}/{full_name}/aligned/{full_name}.trimmed.bam"


            # Map to reference
            echo "mapping ..."
            bwa mem -t 4 {config['reference']} {t_map.inputs[0]} {t_map.inputs[1]}  \
            | samtools view -F 4 -Sb -@ 4 \
            | samtools sort -@ 4 -T {full_name}.align -o $mapped


            # Rename region ids
            echo "renaming ..."
            samtools addreplacerg -@ 4 -r "ID:{full_name}" -o $renamed $mapped
            rm $mapped


            # Trim primers and overall quality
            echo "trimming ..."
            ivar trim -e -i $renamed -b {config['bed_file']} -p $tmptrimmed -q 30
            rm $renamed 
            rm ${{renamed}}.bai

            # Finally sort
            samtools sort -T {full_name}.trim -o {t_map.outputs['bam']} $tmptrimmed
            samtools index {t_map.outputs['bam']}
            rm $tmptrimmed

            
            {default_end}
            """


        # Consensus
        # TODO: Parametrize
        t_consensus = gwf.target(f"cons_{full_name_clean}",
            inputs = t_map.outputs['bam'],
            outputs = f"{output_base}/{full_name}/consensus/{full_name}.fa",
            memory = '16g',
            walltime = '03:00:00') << \
                f"""

                {conda(config['conda_env'])}

                mkdir -p {output_base}/{full_name}/consensus

                samtools mpileup -A -Q 0 -d 0 {output_base}/{full_name}/aligned/{full_name}.sorted.trimmed.bam | ivar consensus -q 30 -t 0.8 -p {output_base}/{full_name}/consensus/{full_name}.fa -m 10 -n N


                {default_end}
                """


        # Pangolin 
        t_pangolin = gwf.target(f"pang_{full_name_clean}",
            inputs = [t_consensus.outputs],
            outputs = [f"{output_base}/{full_name}/pangolin",
                       f"{output_base}/{full_name}/pangolin/{full_name}_pangolin.csv"])
        t_pangolin << \
            f"""

            mkdir -p {t_pangolin.outputs[0]}


            singularity run {config['singularity_images']}/pangolin_latest.sif \
                pangolin {t_pangolin.inputs[0]} \
                    --outdir {t_pangolin.outputs[0]}






            # Prefix header row with #, and end with header for full_name
            cat {t_pangolin.outputs[0]}/lineage_report.csv \
            | head -n 1 \
            | awk '{{ print "#" $0 ",full_name" }}' \
            > {t_pangolin.outputs[1]}

            # End result row with full_name
            cat {t_pangolin.outputs[0]}/lineage_report.csv \
            | tail -n 1 \
            | awk -v sam={full_name} '{{ print $0 "," sam }}' \
            >> {t_pangolin.outputs[1]}


            rm {t_pangolin.outputs[0]}/lineage_report.csv


            {default_end}

            """


        # Nextclade
        t_nextclade = gwf.target(f"next_{full_name_clean}",
            inputs = [t_consensus.outputs],
            outputs = {'dir': f"{output_base}/{full_name}/nextclade",
                       'tab': f"{output_base}/{full_name}/nextclade/{full_name}_nextclade.tab"},
            memory = '4g') 
        t_nextclade << \
            f"""

            mkdir -p {t_nextclade.outputs['dir']}

            # should always print header
            singularity run {config['singularity_images']}/nextclade_latest.sif \
                    nextclade.js \
                        --input-fasta {t_nextclade.inputs[0]} \
                        --output-tsv {t_nextclade.outputs['tab']}.tmp

            ./scripts/dos2unix {t_nextclade.outputs['tab']}.tmp



            echo "catting ..."
            # Put full name in front of all columns, and remove the header row.
            cat {t_nextclade.outputs['tab']}.tmp \
            | awk -v sam={full_name} '{{ print sam "\\t" $0 }}' \
            | grep -vP "seqName\\tclade\\tqc.overallScore" \
            > {t_nextclade.outputs['tab']} || echo -e "{full_name}" > {t_nextclade.outputs['tab']}
            # Forcing graceful exit, because an empty file exits one.


            rm {t_nextclade.outputs['tab']}.tmp

            {default_end}

            """




        #break # Run a single sample only (for testing)

        




print()
