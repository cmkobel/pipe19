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
~~~~~~~~~~~~~~~~~|_|~~~|_|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

""")

"""



"""


input_list_file = "input_list.tab"
batches_done_file = "other/batches_done.tab"
input_base = "input"


df = pd.read_table(input_list_file, sep="\t", names = ["batch", "path", "method"])
batches_done = pd.read_table(batches_done_file, sep = "\t", names = ["batch"])




print(df)
print()
print("//")



for index, row in df.iterrows():

    # if row["bath"] in batches_done:
        #continue

    if row["method"] == "formatA":
        mads_year = 20 
    else:
        mads_year = 21

    print(row['batch'], row['path'])

    try:
        command = f"singularity run docker://rocker/tidyverse Rscript scripts/parse_path.r {row['batch']} {mads_year} {row['path']} TRUE > {input_base}/{row['batch']}.tab" 
        subprocess.run(command, shell = True, check = True)
    except subprocess.CalledProcessError as e:
        print(f"\nAn error occured while initializing {row['batch']}:\n", e)
        sys.exit()


