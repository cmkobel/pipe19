from gwf import Workflow

import pandas as pd
import glob

gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
    "account": "clinicalmicrobio"
})


#input_list_file = "../input_list.tab"
input_glob = "../input/*.tab"
output_base = "output"

# TODO: Use the input files in ../input/* instead of the input_list.tab file. Using the latter means that the batch pipelin will run all over each time the input_list.tab file is touched.
#df = pd.read_table(input_list_file, sep="\t", names = ["batch", "path", "method"], comment = "#")

df = glob.glob(input_glob)


print(df)
print("//")
print()

batch_done_list = []


# Iterate over each line in input_list_file
#for index, row in df.iterrows():
for index, input_file in enumerate(df):
    row = pd.read_table(input_file, sep = "\t")
    print()

    #prefix = f"{row['batch'][0]}"
    prefix = input_file.split("/")[-1].split(".")[0]
    print(f"Batch {index}: {prefix}")
    # if row["bath"] in blacklist:
    #    continue
    print()


    target0 = gwf.target(f"b0_integ_init_{prefix}",
        inputs = input_file,
        outputs = [f"{output_base}/{prefix}/{prefix}_input.tab",
                   f"{output_base}/{prefix}/{prefix}_nextclade.tab",
                   f"{output_base}/{prefix}/{prefix}_pangolin.csv",
                   f"{output_base}/{prefix}/{prefix}_integrated_init.tsv"])
    target0 << \
        f"""

        mkdir -p {output_base}/{prefix}

        echo "input ..."
        cp ../input/{prefix}.tab {target0.outputs[0]}


        echo "nextclade ..."
        cat ../output/{prefix}.*/nextclade/*.tab > {target0.outputs[1]}


        echo "pangolin ..."
        cat ../output/{prefix}.*/pangolin/*.csv > {target0.outputs[2]}


        # Integrate input-pangolin-nextclade files, before joining patient-data
        singularity run ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/integrate_batch_init.r {prefix} {target0.outputs[0]} {target0.outputs[1]} {target0.outputs[2]} {target0.outputs[3]}
            # Rscript args:                               1                    2                    3                    4                    5


        # Cat all integrated_init together
        # TODO: Make it as a regular target with a real Rscript.
        cat output/*/*integrated_init.tsv > integrated_init.tsv

        singularity run ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/fix_tsv.r integrated_init.tsv




        """


    target1 = gwf.target(f"b1_integ_pt___{prefix}",
        inputs = target0.outputs,
        outputs = [f"{output_base}/{prefix}/{prefix}_integrated.tsv",
                   f"{output_base}/{prefix}/{prefix}_sample_sheet.tsv",
                   f"{output_base}/{prefix}/{prefix}_cp_consensus.sh",
                   f"{output_base}/{prefix}/{prefix}_cp_raw.sh"]) #f"{output_base}/{prefix}/{prefix}_upload.tar.gz"]
    target1 << \
        f"""



        



        # TODO: put the catting and compression into another script
        # This will make the surveillance of dangerous variants faster.


        # This file joins everything together, compresses and the files and produces a metadata file for upload
        singularity run ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/integrate_batch.r {prefix} {target1.inputs[3]} "mads/latest/*.csv" {target1.outputs[0]} {target1.outputs[1]}
        # Rscript args:                             1                   2                   3                    4                    5


        # Cat all integrated together
        # TODO: Make it as a regular target with a real Rscript
        cat output/*/*integrated.tsv > integrated.tsv

        singularity run ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/fix_tsv.r integrated.tsv




        """

    # Add the output files for the last per-batch target
    batch_done_list += target1.outputs



    # Copy, compress, and later: upload
    target2 = gwf.target(f"b2_gzip_______{prefix}",
        inputs = target1.outputs,
        outputs = [f"{output_base}/{prefix}/{prefix}_fasta_upload.tar.gz",
                   f"{output_base}/{prefix}/{prefix}_raw_upload.tar.gz",
                   f"{output_base}/{prefix}/{prefix}_compression_done.flag"],
        memory = '4g',
        walltime = '04:00:00')
    target2 << \
        f"""

        #mkdir -p {output_base}/{prefix}/compress
        

        mkdir -p {output_base}/{prefix}/consensus_copy
        mkdir -p {output_base}/{prefix}/raw_copy

        # copy consensus
        echo "copying consensus ..."
        bash {target2.inputs[2]}

        # copy raw fastqs
        echo "copying raw ..."
        bash {target2.inputs[3]}


        # compress
        echo "compressing fasta..."
        cd {output_base}/{prefix}/consensus_copy
        tar -czvf ../../../{target2.outputs[0]} *



        echo "compressing raw..."
        cd ../raw_copy
        tar -czvf ../../../{target2.outputs[1]} *



        echo "touching final flag ..."
        touch ../../../{target2.outputs[2]}



        # TODO: consider removing the raw and consensus copies


        """











    # break # Debug for testing a single batch



#print("batch_done_list", batch_done_list)


target3 = gwf.target(f"b2_coll_all",
    inputs = batch_done_list,
    outputs = "all_batches_integrated.tsv")
target3 << \
    f"""
    singularity run ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/collect_batches.r {target2.outputs}

    """

