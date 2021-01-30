# pipe19

Large scale sequencing of SARS-CoV-2 in Region Midtjylland, Denmark


## Pipeline overview
We use the Artic V3 primers to amplify 98 overlapping amplicons. The amplicons are then fragmented and we sequence each isolate with QiaSeq on Illumina NextSeq and NovaSeq platforms at the local Department of Molecular Medicine.  
This pipeline accepts the demultiplexed data from this sequencing protocol and generally follows the iVar consensus pipeline:

**Steps:**
* Trim index-adaptors fromreads with **trim-galore**.
* Map reads with **bwa mem**.
* Trim primer regions from alignment with **iVar trim**.
* Call consensus with **iVar consensus** (read-depth >= 10, base-frequency >= 0.8).
* Call lineage with **Pangolin**.
* Call clade with **Nextclade**.

For each batch of samples, the following is done:
* Integration of file metadata, pangolin and nextclade calls.
* Compression of raw- and consensus data, tailored for upload.
* Generation of a QC/surveillance report.

**Dependencies:**
* Singularity
* Python 3.7
  * yaml (pyyaml)
  * pandas
  * gwf (gwf-org)

For a complete list of dependencies (and their subdependants), see the conda_env.yml.
  
**Setup:**
Use the config.yaml for setting the input/output directories and the like.
make sure you install all dependencies listed in conda_env.yml. E.g. with the following command: `conda env create -f conda_env.yml`
  
