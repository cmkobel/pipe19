from gwf import Workflow

import pandas as pd
import glob
import yaml

gwf = Workflow(defaults={
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
    "account": "clinicalmicrobio"
})


#input_list_file = "../input_list.tab"
input_glob = "../input/*.tab"
output_base = "output"


with open('../config.yml') as file:
    config = yaml.full_load(file)

#df = pd.read_table(input_list_file, sep="\t", names = ["batch", "path", "method"], comment = "#")

# df = sorted(glob.glob(input_glob))
# 
# print(df)
# print("//")
# print()


input_list = pd.read_table("../" + config['input_list_file'], sep="\t", names = ["batch", "path", "format_specifier"], comment = "#")
input_list["input_file"] = ("../" + config["input_base"] + "/" + input_list["batch"] + ".tab")
print(input_list)
print("//")
print()


batch_done_list = []

#batch_list = []

#for index, row in df.iterrows():
for index, input_list_row in input_list.iterrows():
#for index, input_file in enumerate(df):
    prefix = input_list_row["batch"]
    input_file = input_list_row["input_file"]

    input_table = pd.read_table(input_file, sep = "\t")


    # A little bit of input validation.
    il_batch = str(input_list_row["batch"])
    it_batch = str(set(input_table["batch"]))
    if il_batch not in it_batch:
        raise Exception(f"Discrepancy between content of input list batch ({il_batch}) and input table ({it_batch})")


    print()

    #prefix = f"{row['batch'][0]}"
    #prefix = input_file.split("/")[-1].split(".")[0]
    #batch_list.append(prefix)
    print(f"Batch {index}: {prefix} ({input_file})")
    # if row["batch"] in blacklist:
    #    continue


    target0 = gwf.target(f"init_{prefix}",
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
        singularity run --cleanenv ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/integrate_batch_init.r {prefix} {target0.outputs[0]} {target0.outputs[1]} {target0.outputs[2]} {target0.outputs[3]}
            # Rscript args:                               1                    2                    3                    4                    5


        # Cat all integrated_init together
        # TODO: Make it as a regular target with a real Rscript.
        cat output/*/*integrated_init.tsv > integrated_init.tsv

        singularity run --cleanenv ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/fix_tsv.r integrated_init.tsv




        """

    target0B = gwf.target(f"rprt_{prefix}",
        inputs = target0.outputs,
        outputs = [f"{output_base}/{prefix}/{prefix}_qc_plates_A.pdf",
                   f"{output_base}/{prefix}/{prefix}_qc_plates_B.pdf"])
    target0B << \
        f"""

        singularity run --cleanenv ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/batch_qc.r {target0B.inputs[3]} {prefix} {target0B.outputs[0]} {target0B.outputs[1]} 
            # Rscript args:                               1        2                     3                     4


        """


    target1 = gwf.target(f"pati_{prefix}",
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
        singularity run --cleanenv ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/integrate_batch.r {prefix} {target1.inputs[3]} "mads/latest/*.csv" {target1.outputs[0]} {target1.outputs[1]}
        # Rscript args:                             1                   2                   3                    4                    5


        # Cat all integrated together
        # TODO: Make it as a regular target with a real Rscript
        cat output/*/*integrated.tsv > integrated.tsv

        singularity run --cleanenv ~/faststorage/singularity_images/tidyverse_latest.sif \
            Rscript scripts/fix_tsv.r integrated.tsv




        """

    # Add the output files for the last per-batch target
    batch_done_list += target1.outputs



    # Copy, compress, and later: upload
    target2 = gwf.target(f"gzip_{prefix}",
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

bs = "\n"

#highest_batch_id = sorted(batch_list)[-1]
highest_batch_id = max(input_list["batch"])
print("highest_batch_id", highest_batch_id)
mail_list = open("mail_list.txt", "r").read()

target3 = gwf.target(f"b3_report",
    inputs = batch_done_list,
    outputs = [f"rmarkdown/flags/sent_{highest_batch_id}.flag",
               f"rmarkdown/old_reports/{highest_batch_id}_seq_report.html",],
    memory = '2g')
target3 << \
    f"""

    # Give some slack so a submission can be cancelled
    sleep 60


    # Make sure the old report is cleared if singularity fails without error.
    rm -f rmarkdown/seq_report.html ready.RDS



    # Generate report
    singularity run --cleanenv docker://marcmtk/sarscov2_seq_report \
        render.r rmarkdown/seq_report.Rmd
    


    # {mail_list}
    mail -v -s "Automail: SARS-CoV-2 rapport" -a rmarkdown/seq_report.html {mail_list} <<< "Autogenereret rapport over SARS-CoV-2 i Region Midtjylland (sundhedssporet) til og med sekventeringsbatch-id: {highest_batch_id}

Se vedhæftede html-fil.


__

Klinisk Mikrobioinformatisk Enhed
carkob@rm.dk, marc.nielsen@rm.dk
Klinisk Mikrobiologi ▪ Region Midtjylland
Aarhus Universitetshospital
Palle Juul-Jensens Boulevard 99 ▪ DK-8200 Aarhus

"

    # Simple way of keeping track of whether each new version of the report has been sent.
    mkdir -p rmarkdown/flags
    touch {target3.outputs[0]}


    # Backup the reports
    mkdir -p rmarkdown/old_reports
    cp rmarkdown/seq_report.html {target3.outputs[1]}


    """

target4 = gwf.target(f"b4_voc_list",
    inputs = f"rmarkdown/flags/sent_{highest_batch_id}.flag",
    outputs = f"rmarkdown/old_reports/{highest_batch_id}_voc_list.html",
    memory = '2g')
target4 << f"""

    # Give some slack so a submission can be cancelled
    sleep 60


    # Make sure the old report is cleared if singularity fails without error.
    rm -f rmarkdown/voc_list.html

    # Generate report
    singularity run --cleanenv docker://marcmtk/sarscov2_seq_report \
        render.r rmarkdown/voc_list.Rmd
        
    # Backup the reports
    mkdir -p rmarkdown/old_reports
    cp rmarkdown/voc_list.html {target4.outputs}


    """

mail_list_variant_status = open("mail_list_variant_status.txt", "r").read()
#{mail_list_variant_status}
target5 = gwf.target(f"b5_variant_status",
    inputs = f"rmarkdown/flags/sent_{highest_batch_id}.flag",
    outputs = f"rmarkdown/old_reports/{highest_batch_id}_variant_status.html",
    memory = '2g')
target5 << f"""

    # Give some slack so a submission can be cancelled
     sleep 60


    # Make sure the old report is cleared if singularity fails without error.
    rm -f rmarkdown/variant_status.html

    # Generate report
    singularity run --cleanenv docker://marcmtk/sarscov2_seq_report \
        render.r rmarkdown/variant_status.Rmd
        
    mail -v -s "Automail: SARS-CoV-2 variant status" -a rmarkdown/variant_status.html {mail_list_variant_status} <<< "Autogenereret status over SARS-CoV-2 varianter i Region Midtjylland (sundhedssporet) til og med sekventeringsbatch-id: {highest_batch_id}

Se vedhæftede html-fil.


__

Klinisk Mikrobioinformatisk Enhed
carkob@rm.dk, marc.nielsen@rm.dk
Klinisk Mikrobiologi ▪ Region Midtjylland
Aarhus Universitetshospital
Palle Juul-Jensens Boulevard 99 ▪ DK-8200 Aarhus

"
        
    # Backup the reports
    mkdir -p rmarkdown/old_reports
    cp rmarkdown/variant_status.html {target5.outputs}


    """
