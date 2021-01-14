from gwf import Workflow

from templates import *
import pandas as pd

gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
    "account": "clinicalmicrobio",
    "memory": '2g'
})


input_list_file = "../input_list.tab"
output_base = "output"

df = pd.read_table(input_list_file, sep="\t", names = ["batch", "plate", "path", "method"], comment = "#")


print(df)
print("//")
print()


# Iterate over each line in input_list_file
for index, row in df.iterrows():
    print()
    # if row["bath"] in b

    print(f"Batch {index}: {row['batch']}")

    prefix = f"{row['batch']}.{row['plate']}"


    target0 = gwf.target(f"b0_coll_{prefix}",
        inputs = input_list_file,
        outputs = [f"{output_base}/{prefix}/{prefix}_input.tab",
                   f"{output_base}/{prefix}/{prefix}_nextclade.tab",
                   f"{output_base}/{prefix}/{prefix}_pangolin.csv"])
    target0 << \
        f"""

        mkdir -p {output_base}/{prefix}

        echo "input ..."
        cp ../input/{prefix}.tab {target0.outputs[0]}


        echo "nextclade ..."
        cat ../output/{prefix}.*/nextclade/*.tab > {target0.outputs[1]}


        echo "pangolin ..."
        cat ../output/{prefix}.*/pangolin/*.csv > {target0.outputs[2]}

        """


    target1 = gwf.target(f"b1_R____{prefix}",
        inputs = target0.outputs,
        outputs = f"{output_base}/{prefix}/{prefix}_sample_sheet.tsv", 
        memory = '4g')
    target1 << \
        f"""

        # This file joins everything together, compresses and the files and produces a metadata file for SSI
        singularity run ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/integrate_batch.r

        """




    break

