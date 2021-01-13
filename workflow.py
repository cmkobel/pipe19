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

        print(full_name)

        #print(" ", [batch_row['path'] + i for i in (batch_row['R1'] + batch_row['R2']).split(" ")]); exit()
        target_0cat = gwf.target(f"_0cat_{full_name}".replace(".", "_"),
            inputs = [batch_row['path'] + i for i in batch_row['R1'].split(" ")] +
                     [batch_row['path'] + i for i in batch_row['R2'].split(" ")], 
            outputs = [f"{output_base}/{full_name}/trimmed/{full_name}_val_1.fq.gz",
                       f"{output_base}/{full_name}/trimmed/{full_name}_val_2.fq.gz"],
            cores = 4) << \
                f"""

                {conda("ivar-inpipe")}
                
                mkdir -p {output_base}/{full_name}/trimmed
                

                # Cat

                cat {" ".join([batch_row['path'] + i for i in batch_row['R1'].split(" ")])} > {output_base}/{full_name}/trimmed/{full_name}_R1{batch_row['extension']}
                cat {" ".join([batch_row['path'] + i for i in batch_row['R2'].split(" ")])} > {output_base}/{full_name}/trimmed/{full_name}_R2{batch_row['extension']}


                # Trim

                trim_galore --paired --fastqc --cores 4 --gzip -o {output_base}/{full_name}/trimmed --basename {full_name} {output_base}/{full_name}/trimmed/{full_name}_R1{batch_row['extension']} {output_base}/{full_name}/trimmed/{full_name}_R2{batch_row['extension']}




        

            """

       
    
                # Map 
        target_1map = gwf.target(f"_1map_{full_name}".replace(".", "_"),
            inputs = target_0cat.outputs,
            outputs = [f"{output_base}/{full_name}/aligned/{full_name}.sorted.bam",
                       f"{output_base}/{full_name}/aligned/{full_name}.sorted.trimmed.bam"],
            cores = 4) << \
                f"""

                {conda("ivar-inpipe")} # TODO: remove bwa from pipe19_a

                mkdir -p {output_base}/{full_name}/aligned


                # Map to reference
                echo "mapping ..."
                bwa mem -t 4 {config['reference']} {output_base}/{full_name}/trimmed/{full_name}_val_1.fq.gz {output_base}/{full_name}/trimmed/{full_name}_val_2.fq.gz \
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


        break

