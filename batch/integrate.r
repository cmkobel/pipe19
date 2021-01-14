# First, call collect.sh to collect all the results together

# Then import it

rm(list = ls())
library(tidyverse)
# setwd("~/GenomeDK/ClinicalMicrobio/faststorage/pipe19/batch")


df_input = read_tsv("collected/input.tab") %>% 
    mutate(full_name = paste0(batch, ".", plate, ".", moma_serial, "_", sample_name)) %>% 
    select(full_name, everything())


df_nextclade = read_tsv("collected/nextclade.tab", comment = "#", skip = 1,
                        col_names = c("seqName", "clade", "qc.overallScore", "qc.overallStatus", "totalGaps", "totalInsertions", "totalMissing", "totalMutations", "totalNonACGTNs", "totalPcrPrimerChanges", "substitutions", "deletions", "insertions", "missing", "nonACGTNs", "pcrPrimerChanges", "aaSubstitutions", "totalAminoacidSubstitutions", "aaDeletions", "totalAminoacidDeletions", "alignmentEnd", "alignmentScore", "alignmentStart", "qc.missingData.missingDataThreshold", "qc.missingData.score", "qc.missingData.status", "qc.missingData.totalMissing", "qc.mixedSites.mixedSitesThreshold", "qc.mixedSites.score", "qc.mixedSites.status", "qc.mixedSites.totalMixedSites", "qc.privateMutations.cutoff", "qc.privateMutations.excess", "qc.privateMutations.score", "qc.privateMutations.status", "qc.privateMutations.total", "qc.snpClusters.clusteredSNPs", "qc.snpClusters.score", "qc.snpClusters.status", "qc.snpClusters.totalSNPs", "errors", "full_name"))

df_pangolin = read_delim("collected/pangolin.csv", delim = ",", comment = "#", skip = 1,
                         col_names = c("taxon", "lineage", "probability", "pangoLEARN_version", "status", "note", "full_name"))




integrated = left_join(df_input, df_pangolin, by = "full_name") %>% 
    left_join(df_nextclade, by = "full_name")


integrated %>% write_rds("integrated.rds")


