from gwf import Workflow

import pandas as pd

gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
    "account": "clinicalmicrobio",
    "memory": '2g',
    "walltime": "02:00:00"
})


input_list_file = "../input_list.tab"
output_base = "output"

df = pd.read_table(input_list_file, sep="\t", names = ["batch", "path", "method"], comment = "#")


print(df)
print("//")
print()

batch_done_list = []


# Iterate over each line in input_list_file
for index, row in df.iterrows():
    print()
    # if row["bath"] in b

    print(f"Batch {index}: {row['batch']}")

    prefix = f"{row['batch']}"


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
        outputs = [f"{output_base}/{prefix}/{prefix}_integrated.tsv",
                   f"{output_base}/{prefix}/{prefix}_sample_sheet.tsv",
                   f"{output_base}/{prefix}/{prefix}_upload.tar.gz"],
        memory = '4g',
        walltime = '04:00:00')
    target1 << \
        f"""

        # mkdir -p {output_base}/{prefix}/raw_copy
        # mkdir -p {output_base}/{prefix}/consensus_copy

        mkdir -p {output_base}/{prefix}/compress




        # This file joins everything together, compresses and the files and produces a metadata file for SSI
        singularity run ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/integrate_batch.r {prefix} {target1.inputs[0]} {target1.inputs[1]} {target1.inputs[2]} "mads/latest/*.csv" {target1.outputs[0]} {target1.outputs[1]} {target1.outputs[2]}
        #                                            1                   2                   3                   4                   5                    6                    7                    8


        # TODO: consider removing the raw and consensus copies




        """

    # Add the output files for the last per-batch target
    batch_done_list += target1.outputs




    # break # Debug for testing a single batch



#print("batch_done_list", batch_done_list)


target2 = gwf.target(f"b2_collect_all",
    inputs = batch_done_list,
    outputs = "all_batches_integrated.tsv")
target2 << \
    f"""
    singularity run ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/collect_batches.r {target2.outputs}

    """

