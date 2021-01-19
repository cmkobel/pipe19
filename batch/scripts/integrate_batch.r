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
file_mads = args[5]
file_integrated_out = args[6]
file_sample_sheet_out = args[7]
file_targz_out = args[8]



if (development_mode) {
    rm(list = ls())
    
    batch = "210108"
    
    # Inputs
    file_input = "output/210108/210108_input.tab"
    file_nextclade = "output/210108/210108_nextclade.tab"
    file_pangolin = "output/210108/210108_pangolin.csv"
    file_mads = "mads/latest/*.csv"
    
    # Outputs
    file_integrated_out = "output/210108/210108_integrated.tsv"
    file_sample_sheet_out = "output/210108/210108_samplesheet.tsv"
    file_targz_out = "output/210108/210108_upload.tar.gz"
}

write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

write("importing input list ...", stderr())
df_input = read_tsv(file_input) %>% 
    mutate(full_name = paste0(batch, ".", plate, ".", moma_serial, "_", sample_name)) %>% 
    select(full_name, everything())

write("importing nextclade ...", stderr())
df_nextclade = read_tsv(file_nextclade, comment = "#", skip = 1,
                        col_names = c("seqName", "clade", "qc.overallScore", "qc.overallStatus", "totalGaps", "totalInsertions", "totalMissing", "totalMutations", "totalNonACGTNs", "totalPcrPrimerChanges", "substitutions", "deletions", "insertions", "missing", "nonACGTNs", "pcrPrimerChanges", "aaSubstitutions", "totalAminoacidSubstitutions", "aaDeletions", "totalAminoacidDeletions", "alignmentEnd", "alignmentScore", "alignmentStart", "qc.missingData.missingDataThreshold", "qc.missingData.score", "qc.missingData.status", "qc.missingData.totalMissing", "qc.mixedSites.mixedSitesThreshold", "qc.mixedSites.score", "qc.mixedSites.status", "qc.mixedSites.totalMixedSites", "qc.privateMutations.cutoff", "qc.privateMutations.excess", "qc.privateMutations.score", "qc.privateMutations.status", "qc.privateMutations.total", "qc.snpClusters.clusteredSNPs", "qc.snpClusters.score", "qc.snpClusters.status", "qc.snpClusters.totalSNPs", "errors", "full_name"))


df_pangolin = read_delim(file_pangolin, delim = ",", comment = "#", skip = 1,
                         col_names = c("taxon", "lineage", "probability", "pangoLEARN_version", "status", "note", "full_name"))


write("importing df_mads ...", stderr())
df_mads = read_csv(Sys.glob(file_mads) %>% sort %>% tail(1), locale = locale(encoding = "WINDOWS-1252")) %>% 
#names(df_mads[5]) = "proevenr"
#write(paste("these are the names of df_mads", names(df_mads)), stderr())
#df_mads = df_mads %>% 
    rename(koen = `K<U+00D8>N`, # Hej Marc. Ved ved ikke hvorfor der er problemer med indkodningen af kolonne-navne. Derfor laver jeg manuelle oversættelser her.
           proevenr = `pr<U+00F8>venr`,
           proevekategori = `pr<U+00F8>vekategori`) %>% 
    mutate_at(vars(afsendt, modtaget), lubridate::dmy) %>%
    ##mutate_at(vars(klokken), lubridate::hm) %>% 
    mutate(ct = coalesce(X16, X19, X22, X25, X28, X31))

write(paste("these are the names of df_mads", paste(names(df_mads), collapse = ", ")), stderr())
    

# Join everything together
df_integrated = left_join(df_input, df_pangolin, by = "full_name") %>% 
    left_join(df_nextclade, by = "full_name") %>% 
    left_join(df_mads, by = c("sample_name" = "proevenr"))


# Save the full table
write(paste("writing df_integrated to", file_integrated_out), stderr())
df_integrated %>% write_tsv(file_integrated_out)



# Make the file ready for export to the government
## Filter for samples only, and select/rename the columns of interest.
# ”sample_id;cpr;sampling_date;kma_id;raw_filename;consensus_filename”.
write(paste("writing sample sheet to", file_sample_sheet_out), stderr())
df_sample_sheet = df_integrated %>% filter(type == "sample") %>% 
    rowwise() %>% 
    mutate(kma_id = "6620320",
           full_name = paste0(batch, ".", plate, ".", moma_serial, "_", sample_name),
           raw_filename = paste0(full_name, "_R", c(1, 2), ".fastq.gz", collapse = " "),
           consensus_filename = paste0(full_name, ".fa", collapse = " "),
           platform = "illumina qiaseq") %>% 
    ungroup() %>% 
    select(sample_id = sample_name, full_name, cpr = `cprnr.`, sampling_date = afsendt, kma_id, raw_filename, consensus_filename, platform, ct) # Consider including Ydernr/SKSnr

df_sample_sheet %>%  
    select(-full_name) %>% 
    write_tsv(file_sample_sheet_out)
  



# tar the files together,
# It makes sense to do it from here, because the relevant metadata is already loaded in the environment.
target_copy = paste0("output/", batch, "/compress")
#target_raw = paste0("output/", batch, "/raw_copy/")
#target_consensus = paste0("output/", batch, "/consensus_copy/")

  

# Generate a table that has the commands. Then paste them to a system-command

# Copy raw data
write("commanding raws ...", stderr())
command_raw = df_sample_sheet %>% 
    select(full_name, sample_id, raw_filename) %>% 
    rowwise() %>% 
    mutate(raw_filename_splitted = str_split(raw_filename, " "),
           source_files = paste0("../output/", full_name, "/trimmed_reads/", raw_filename_splitted, collapse = " "),
           target_dir = target_copy) %>% 
    ungroup() %>% 
    
    transmute(command = paste("cp", source_files, target_dir))

command_raw %>% write_tsv(paste0("output/", batch, "/cp_raw.log"))

if (!development_mode) {
  system(command_raw %>% pull(command) %>% paste(collapse = "; "))
}


# Copy Consensus data
write("commanding consensus' ...", stderr())
command_consensus = df_sample_sheet %>% 
    select(full_name, sample_id, consensus_filename) %>% 
    rowwise() %>% 
    mutate(source_file = paste0("../output/", full_name, "/consensus/", consensus_filename),
         target_dir = target_copy) %>% 
    ungroup() %>% 
  
    transmute(command = paste("cp", source_file, target_dir))

command_consensus %>% write_tsv(paste0("output/", batch, "/cp_consensus.log"))

if (!development_mode) {
    system(command_consensus %>% pull(command) %>% paste(collapse = "; "))
}

# Tar the files
# Could also be called from the workflow. But here it is possibly easier.
# TODO: Consider deleting some intermediary files

write("compressing ...", stderr())
system(paste0("cd ", target_copy, "; tar -czvf ../../../", file_targz_out, " *"))







       
# upload the files via FTP

write("\nDone for now.", stderr())




