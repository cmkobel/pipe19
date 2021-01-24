# pipe19

Large scale sequencing of SARS-CoV-2 in Region Midtjylland, Denmark


## Pipeline overview
Swabs are purified, and the RNA is transcribed to cDNA. We use the Artic V3 primers to amplify 98 overlapping amplicons. The amplicons are then fragmented and we sequence each isolate with QiaSeq on Illumina NextSeq and NovaSeq platforms at the local Department of Molecular Medicine.  
This pipeline accepts the demultiplexed data from this sequencing protocol and generally follows the iVar consensus pipeline:

**Steps:**
* Trim index-adaptors fromreads with ***trim-galore***.
* Map reads with ***bwa mem***.
* Trim primer regions from alignment with ***iVar trim***.
* Call consensus with ***iVar consensus*** (read-depth >= 10, base-frequency >= 0.8).
* Call lineage with ***Pangolin***.
* Call clade with ***Nextclade***.
