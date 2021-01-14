# First, call collect.sh to collect all the results together

# Then import it


library(tidyverse)
# setwd("~/GenomeDK/ClinicalMicrobio/faststorage/pipe19/batch")



args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())


batch = args[1]
file_input = args[2]
file_nextclade = args[3]
file_pangolin = args[4]


devel = T

if (devel) {
    rm(list = ls())
    batch = "210108.3471"
    file_input = "output/210108.3471/210108.3471_input.tab"
    file_nextclade = "output/210108.3471/210108.3471_nextclade.tab"
    file_pangolin = "output/210108.3471/210108.3471_pangolin.csv"
    file_mads = "mads/latest/*.csv"
    file_out = "output/210108.3471/210108.3471_samplesheet.tsv"
}


df_input = read_tsv(file_input) %>% 
    mutate(full_name = paste0(batch, ".", plate, ".", moma_serial, "_", sample_name)) %>% 
    select(full_name, everything())


df_nextclade = read_tsv(file_nextclade, comment = "#", skip = 1,
                        col_names = c("seqName", "clade", "qc.overallScore", "qc.overallStatus", "totalGaps", "totalInsertions", "totalMissing", "totalMutations", "totalNonACGTNs", "totalPcrPrimerChanges", "substitutions", "deletions", "insertions", "missing", "nonACGTNs", "pcrPrimerChanges", "aaSubstitutions", "totalAminoacidSubstitutions", "aaDeletions", "totalAminoacidDeletions", "alignmentEnd", "alignmentScore", "alignmentStart", "qc.missingData.missingDataThreshold", "qc.missingData.score", "qc.missingData.status", "qc.missingData.totalMissing", "qc.mixedSites.mixedSitesThreshold", "qc.mixedSites.score", "qc.mixedSites.status", "qc.mixedSites.totalMixedSites", "qc.privateMutations.cutoff", "qc.privateMutations.excess", "qc.privateMutations.score", "qc.privateMutations.status", "qc.privateMutations.total", "qc.snpClusters.clusteredSNPs", "qc.snpClusters.score", "qc.snpClusters.status", "qc.snpClusters.totalSNPs", "errors", "full_name"))


df_pangolin = read_delim(file_pangolin, delim = ",", comment = "#", skip = 1,
                         col_names = c("taxon", "lineage", "probability", "pangoLEARN_version", "status", "note", "full_name"))


df_mads = read_csv(Sys.glob(file_mads) %>% sort %>% tail(1), locale = locale(encoding = "WINDOWS-1252")) %>% 
    mutate_at(vars(afsendt, modtaget), lubridate::dmy) %>%
    #mutate_at(vars(klokken), lubridate::hm) %>% 
    mutate(ct = coalesce(X16, X19, X22, X25, X28, X31, X34))


# Join everything together
df_integrated = left_join(df_input, df_pangolin, by = "full_name") %>% 
    left_join(df_nextclade, by = "full_name") %>% 
    left_join(df_mads, by = c("sample_name" = "prøvenr"))


# Make the file ready for export and 




# Collect and compress files




