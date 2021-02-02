# First, call collect.sh to collect all the results together

# Then import it


library(tidyverse)
development_mode = F
# setwd("~/GenomeDK/ClinicalMicrobio/faststorage/pipe19/batch"); development_mode = T



args = commandArgs(trailingOnly=TRUE)



batch = args[1]
file_input = args[2]
file_nextclade = args[3]
file_pangolin = args[4]
file_integrated_out = args[5]



if (development_mode) {
    rm(list = ls())
    
    batch = "210108"
    
    # Inputs
    file_input = "output/210108/210108_input.tab"
    file_nextclade = "output/210108/210108_nextclade.tab"
    file_pangolin = "output/210108/210108_pangolin.csv"

    # Outputs
    file_integrated_out = "output/210108/210108_integrated_init.tsv"
}

write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

write("importing input list ...", stderr())
raw_input = read_tsv(file_input)

df_input = raw_input %>% 
    
    mutate(sample_name_prefix = str_sub(raw_sample_name, 1, 2), # first two characters
           sample_name_suffix = str_sub(raw_sample_name, -6),   # last six characters
           raw_full_name = paste0(batch, ".", plate, ".", moma_serial, "_", raw_sample_name),
           sample_name_prefix_converted = recode(sample_name_prefix,
                                                 "87" = paste0("L"),
                                                 "88" = paste0("R"),
                                                 "89" = paste0("I"),
                                                 "90" = paste0("V"),
                                                 "96" = paste0("P"))) %>% 
          
    rowwise() %>% 
    mutate(ya_sample_name = if_else(type == "sample", # "year_agnostic_sample_name"
                                         paste0(sample_name_prefix_converted, sample_name_suffix),
                                         paste0(raw_sample_name))) %>% 
    ungroup() %>% 


    select(raw_full_name, batch, plate, moma_serial, raw_sample_name, ya_sample_name, type, extension, source_project, path, R1, R2)







write("importing nextclade ...", stderr())
df_nextclade = read_tsv(file_nextclade, comment = "#", skip = 1,
                        col_names = c("seqName", "clade", "qc.overallScore", "qc.overallStatus", "totalGaps", "totalInsertions", "totalMissing", "totalMutations", "totalNonACGTNs", "totalPcrPrimerChanges", "substitutions", "deletions", "insertions", "missing", "nonACGTNs", "pcrPrimerChanges", "aaSubstitutions", "totalAminoacidSubstitutions", "aaDeletions", "totalAminoacidDeletions", "alignmentEnd", "alignmentScore", "alignmentStart", "qc.missingData.missingDataThreshold", "qc.missingData.score", "qc.missingData.status", "qc.missingData.totalMissing", "qc.mixedSites.mixedSitesThreshold", "qc.mixedSites.score", "qc.mixedSites.status", "qc.mixedSites.totalMixedSites", "qc.privateMutations.cutoff", "qc.privateMutations.excess", "qc.privateMutations.score", "qc.privateMutations.status", "qc.privateMutations.total", "qc.snpClusters.clusteredSNPs", "qc.snpClusters.score", "qc.snpClusters.status", "qc.snpClusters.totalSNPs", "errors", "raw_full_name"))


df_pangolin = read_delim(file_pangolin, delim = ",", comment = "#", skip = 1,
                         col_names = c("taxon", "lineage", "probability", "pangoLEARN_version", "status", "note", "raw_full_name"))




# Join everything together
df_integrated = left_join(df_input, df_pangolin, by = "raw_full_name") %>% 
    left_join(df_nextclade, by = "raw_full_name")


# Save the full table
write(paste("writing df_integrated to", file_integrated_out), stderr())
df_integrated %>% write_tsv(file_integrated_out)





