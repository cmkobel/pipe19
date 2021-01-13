#! /usr/bin/env python3
from gwf import *
import glob
import os 
#from cutils import *
import json

import pandas as pd
import subprocess, sys



gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
    "account": "clinicalmicrobio"
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


input_list_file = "input_list.tab"
batches_done_file = "other/batches_done.tab"
input_base = "input"
output_base = "output"


config = {'reference': "../artic/ivar/artic-ncov2019_gh_clone/artic-ncov2019/primer_schemes/SARS-CoV-2/V3/nCoV-2019.reference.fasta",
          'bed_file': "../artic/ivar/artic-ncov2019_gh_clone/artic-ncov2019/primer_schemes/SARS-CoV-2/V3/nCoV-2019.primer.bed"}


df = pd.read_table(input_list_file, sep="\t", names = ["batch", "plate", "path", "method"])
batches_done = pd.read_table(batches_done_file, sep = "\t", names = ["batch"])


def conda(env):
    return f"source ~/miniconda3/etc/profile.d/conda.sh; conda activate {env}"


print("input list file:")
print(df)
print("//")
print()


# Iterate over each line in input_list_file
for index, row in df.iterrows():

    # if row["bath"] in batches_done:
        #continue

    if row["method"] == "formatA":
        mads_year = 20 
    else:
        mads_year = 21

    batch_long = f"{row['batch']}.{row['plate']}" # e.g. 210108.3471

    print(f"Batch {index}: {batch_long}, {row['path']}")

    batch_input_file = f"{input_base}/{batch_long}.tab"

    # Check if the batch_input_file has already been written
    if os.path.exists(batch_input_file):
        print(f" {batch_input_file} has already been written")
    else:
        print(f" creating the file {batch_input_file}")
        try:
            command = f"singularity run docker://rocker/tidyverse Rscript scripts/parse_path.r {row['batch']} {row['plate']} {row['path']} {mads_year} TRUE > {batch_input_file}" 
            subprocess.run(command, shell = True, check = True)
        except subprocess.CalledProcessError as e:
            print(f"\nAn error occured while initializing {row['batch']}:\n", e)
            sys.exit()


    # Now the batch_input_file file has surely been written, and we can start the actual pipeline
    batch_df = pd.read_table(batch_input_file, sep = "\t")

    print(batch_df)
    print("//")

    # TODO: Check that the batch given in input list is similar to the batch specified in the batch_input_file

    for batch_index, batch_row in batch_df.iterrows():

        sample_name = f"{batch_row['sample_name']}"
        full_name = f"{batch_long}.{batch_row['moma_serial']}_{batch_row['sample_name']}"
        full_name_clean = full_name.replace(".", "_")

        print(full_name)

        #print(" ", [batch_row['path'] + i for i in (batch_row['R1'] + batch_row['R2']).split(" ")]); exit()
        t_cat = gwf.target(f"_cat__{full_name_clean}",
            inputs = [batch_row['path'] + i for i in batch_row['R1'].split(" ")] +
                     [batch_row['path'] + i for i in batch_row['R2'].split(" ")], 
            outputs = {'dir': f"{output_base}/{full_name}/trimmed_reads/",
                       'files': [f"{output_base}/{full_name}/trimmed_reads/{full_name}_val_1.fq.gz",
                                 f"{output_base}/{full_name}/trimmed_reads/{full_name}_val_2.fq.gz"]},
            cores = 4)
        t_cat << \
            f"""

            {conda("ivar-inpipe")}
            
            mkdir -p {t_cat.outputs['dir']}
            
            tmp_forward={output_base}/{full_name}/trimmed_reads/{full_name}_R1{batch_row['extension']}
            tmp_reverse={output_base}/{full_name}/trimmed_reads/{full_name}_R2{batch_row['extension']}


            # Cat the reads together
            cat {" ".join(t_cat.inputs[0])} > $tmp_forward
            cat {" ".join(t_cat.inputs[1])} > $tmp_reverse


            # Trim the reads
            trim_galore --paired --fastqc --cores 4 --gzip -o {t_cat.outputs['dir']} --basename {full_name} $tmp_forward $tmp_reverse


            # TODO: Consider removing the catted reads (tmp).



        

            """

       
    
        # Map 
        t_map = gwf.target(f"_map__{full_name_clean}",
            inputs = t_cat.outputs['files'],
            outputs = [f"{output_base}/{full_name}/aligned/{full_name}.sorted.trimmed.bam"],
            cores = 4) << \
                f"""

                {conda("ivar-inpipe")} # TODO: remove bwa from pipe19_a

                mkdir -p {output_base}/{full_name}/aligned


                # Map to reference
                echo "mapping ..."
                bwa mem -t 4 {config['reference']} {output_base}/{full_name}/trimmed_reads/{full_name}_val_1.fq.gz {output_base}/{full_name}/trimmed_reads/{full_name}_val_2.fq.gz \
                | samtools view -F 4 -Sb -@ 4 \
                | samtools sort -@ 4 -T {full_name}.align -o {output_base}/{full_name}/aligned/{full_name}.sorted.tmp.bam


                # Rename region ids
                echo "renaming ..."
                samtools addreplacerg -@ 4 -r "ID:{full_name}" -o {output_base}/{full_name}/aligned/{full_name}.sorted.bam {output_base}/{full_name}/aligned/{full_name}.sorted.tmp.bam
                rm {output_base}/{full_name}/aligned/{full_name}.sorted.tmp.bam


                # Trim alignment
                echo "trimming ..."
                ivar trim -e -i {output_base}/{full_name}/aligned/{full_name}.sorted.bam -b {config['bed_file']} -p {output_base}/{full_name}/aligned/{full_name}.trimmed.bam -q 30
                samtools sort -T {full_name}.trim -o {output_base}/{full_name}/aligned/{full_name}.sorted.trimmed.bam {output_base}/{full_name}/aligned/{full_name}.trimmed.bam
                rm {output_base}/{full_name}/aligned/{full_name}.trimmed.bam



                """


        # Consensus
        t_consensus = gwf.target(f"_cons_{full_name_clean}",
            inputs = t_map.outputs,
            outputs = f"{output_base}/{full_name}/consensus/{full_name}.fa") << \
                f"""

                {conda("ivar-inpipe")}

                mkdir -p {output_base}/{full_name}/consensus

                samtools mpileup -A -Q 0 -d 0 {output_base}/{full_name}/aligned/{full_name}.sorted.trimmed.bam | ivar consensus -q 30 -p {output_base}/{full_name}/consensus/{full_name}.fa -m 10 -n N


                """

        # Pangolin 
        t_pangolin = gwf.target(f"_pang_{full_name_clean}",
            inputs = t_consensus.outputs,
            outputs = [f"{output_base}/{full_name}/pangolin",
                       f"{output_base}/{full_name}/pangolin/{full_name}_pangolin.csv"])
        t_pangolin << \
            f"""

            mkdir -p {t_pangolin.outputs[0]}


            singularity run docker://staphb/pangolin \
                pangolin {t_pangolin.inputs} \
                    --threads 1 \
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


            """



        break

        




