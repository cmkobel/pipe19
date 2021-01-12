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


df = pd.read_table(input_list_file, sep="\t", names = ["batch", "plate", "path", "method"])
batches_done = pd.read_table(batches_done_file, sep = "\t", names = ["batch"])



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
        gwf.target(f"_0cat_{full_name}".replace(".", "_"),
            inputs = [batch_row['path'] + i for i in batch_row['R1'].split(" ")] +
                     [batch_row['path'] + i for i in batch_row['R2'].split(" ")], 
            outputs = [f"{output_base}/{full_name}/catted/{full_name}_R1.{batch_row['extension']}",
                       f"{output_base}/{full_name}/catted/{full_name}_R2.{batch_row['extension']}"]) << \
            f"""

            mkdir -p {output_base}/{full_name}/catted/
            cat {" ".join([batch_row['path'] + i for i in batch_row['R1'].split(" ")])} > {output_base}/{full_name}/catted/{full_name}_R1{batch_row['extension']}
            cat {" ".join([batch_row['path'] + i for i in batch_row['R2'].split(" ")])} > {output_base}/{full_name}/catted/{full_name}_R2{batch_row['extension']}
            

            """

        # Trim 


        # Map 


        # Call consensus


        # Nextclade


        # Pangolin




