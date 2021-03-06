# First, collect all the results together.

# Then import it


library(tidyverse)
development_mode = F
# development_mode = T



args = commandArgs(trailingOnly=TRUE)



batch = args[1]
#file_input = args[2]
#file_nextclade = args[3]
#file_pangolin = args[4]
file_integrated_init = args[2]
file_mads = args[3]
file_integrated_out = args[4]
file_sample_sheet_out = args[5]



if (development_mode) {
  
    development_mode = T
    
    setwd("~/GenomeDK/clinmicrocore/pipe19/batch")
    
    batch = "210223"
    
    # Inputs
    file_integrated_init = "output/210223/210223_integrated_init.tsv"
    file_mads = "mads/latest/*.csv"
    
    # Outputs
    file_integrated_out = "output/210223/210223_integrated.tsv"
    file_sample_sheet_out = "output/210223/210223_sample_sheet.tsv"
}

write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

write("importing integrated_init ...", stderr())
df_integrated_init = read_tsv(file_integrated_init)

# df_input = read_tsv(file_input) %>% 
#     mutate(full_name = paste0(batch, ".", plate, ".", moma_serial, "_", sample_name)) %>% 
#     select(full_name, everything())
# 
# write("importing nextclade ...", stderr())
# df_nextclade = read_tsv(file_nextclade, comment = "#", skip = 1,
#                         col_names = c("seqName", "clade", "qc.overallScore", "qc.overallStatus", "totalGaps", "totalInsertions", "totalMissing", "totalMutations", "totalNonACGTNs", "totalPcrPrimerChanges", "substitutions", "deletions", "insertions", "missing", "nonACGTNs", "pcrPrimerChanges", "aaSubstitutions", "totalAminoacidSubstitutions", "aaDeletions", "totalAminoacidDeletions", "alignmentEnd", "alignmentScore", "alignmentStart", "qc.missingData.missingDataThreshold", "qc.missingData.score", "qc.missingData.status", "qc.missingData.totalMissing", "qc.mixedSites.mixedSitesThreshold", "qc.mixedSites.score", "qc.mixedSites.status", "qc.mixedSites.totalMixedSites", "qc.privateMutations.cutoff", "qc.privateMutations.excess", "qc.privateMutations.score", "qc.privateMutations.status", "qc.privateMutations.total", "qc.snpClusters.clusteredSNPs", "qc.snpClusters.score", "qc.snpClusters.status", "qc.snpClusters.totalSNPs", "errors", "full_name"))
# 
# 
# df_pangolin = read_delim(file_pangolin, delim = ",", comment = "#", skip = 1,
#                          col_names = c("taxon", "lineage", "probability", "pangoLEARN_version", "status", "note", "full_name"))

file_mads_latest = Sys.glob(file_mads) %>% sort %>% tail(1)
write(paste("importing latest df_mads", file_mads_latest, "..."), stderr())
#df_mads = read_csv(file_mads_latest, locale = locale(encoding = "WINDOWS-1252")) %>% 
df_mads = read_csv(file_mads_latest) %>% 
  

    rename(koen = `K<U+00D8>N`, 
           proevenr = `pr<U+00F8>venr`,
           proevekategori = `pr<U+00F8>vekategori`) %>% 
    mutate_at(vars(afsendt, modtaget), lubridate::dmy) %>%
    ##mutate_at(vars(klokken), lubridate::hm) %>% 
    mutate(ct = coalesce(X16, X19, X22, X25, X28, X31, X34, X37)) %>% 
  
    mutate(final_sample_name = proevenr,
           ya_sample_name = paste0(str_sub(final_sample_name, 1, 1), str_sub(final_sample_name, -6)))

write(paste("these are the names of df_mads", paste(names(df_mads), collapse = ", ")), stderr())
    
#df_mads %>% select(X37, ct) %>% View


# emulate conflict with older proevenr
# mads_additional = tibble(ya_sample_name = "R081488",
#                        modtaget = lubridate::ymd("2019-12-30"))
# df_mads = bind_rows(df_mads,
#                     mads_additional)


# Join everything together, taking care of conflicts in the dates..
df_integrated_ranked = left_join(df_integrated_init, df_mads, by = "ya_sample_name") %>% 
  group_by(ya_sample_name) %>% 
  mutate(rank = row_number((desc(modtaget))))
  

# output conflicts for surveillance
df_integrated_ranked %>% 
  filter(!is.na(rank)) %>% 
  filter(length(ya_sample_name) >= 2) %>% 
  ungroup() %>% 
  select(`#batch` = batch, final_sample_name, ya_sample_name, raw_sample_name, modtaget, rank) %>% 
  write_tsv(paste0("output/", batch, "/conflicts.tsv"), col_names = F)
  
# After outputting the date-conflicts, the newest samples can now be filtered for.
df_integrated = df_integrated_ranked %>% 
  filter(rank == 1 | is.na(rank)) %>%  # Always picks the newest sample, when there are conflicts. The ones with NA-value are the ones where no match was found in the mads-extract, possibly because the mads-extract is too old.
  ungroup() %>% 
  select(final_sample_name, everything(), -rank)




# Check that no missing values exist.
missing = df_integrated %>%
  filter(type == "sample") %>% # Comment this to check the behaviour.
  filter(is.na(`cprnr.`))



write(paste("missing cpr:", missing$raw_full_name), stdout())

if (dim(missing)[1] > 0) {
    fail_msg = paste("Missing values in integrated prior to writing for the following raw sample names:\n", paste(missing$raw_full_name, collapse = ", "), "\n\nPlease update to the newest mads extract before continuing.")
    #write(fail_msg, stderr())
    
    fail_on_purpose = T
} else {
    fail_on_purpose = F 
}



#df_integrated_ranked %>% select(modtaget, rank, ya_sample_name) %>% View


# Save the full integrated table
# This file will, later, be concatenated across all batches.
write(paste("writing df_integrated to", file_integrated_out), stderr())
df_integrated %>% write_tsv(file_integrated_out)



# Make the file ready for export to the government
## Filter for samples only, and select/rename the columns of interest.
# ”sample_id;cpr;sampling_date;kma_id;raw_filename;consensus_filename”.
write(paste("writing sample sheet to", file_sample_sheet_out), stderr())
df_sample_sheet = df_integrated %>%
    filter(type == "sample") %>% 
    filter(plate_control_summary != "unsatisfactory") %>%  
    rowwise() %>% 
    mutate(kma_id = "6620320",
           sample_id = raw_sample_name,
           #raw_full_name = paste0(batch, ".", plate, ".", moma_serial, "_", raw_sample_name),
           raw_filename = paste0(raw_full_name, "_R", c(1, 2), ".fastq.gz", collapse = ";"),
           #consensus_filename = paste0(sample_id, ".fa", collapse = " "), # hvorfor er den her collapsed?
            consensus_filename = paste0(sample_id, ".fasta"),
           platform = "illumina") %>% 
    ungroup() %>% 
    select(raw_full_name, sample_id, cpr = `cprnr.`, sampling_date = afsendt, kma_id, raw_filename, consensus_filename, platform, ct) # Consider including Ydernr/SKSnr

df_sample_sheet %>%  
    select(-raw_full_name) %>% 
    write_tsv(file_sample_sheet_out)
  


# tar the files together,
# It makes sense to do it from here, because the relevant metadata is already loaded in the environment.
#target_copy = paste0("output/", batch, "/compress")
target_raw = paste0("output/", batch, "/raw_copy/")
target_consensus = paste0("output/", batch, "/consensus_copy/")

  

# Generate a table that has the commands. Then paste them to a system-command

# Copy raw data
# I have disabled this block because SSI doesn't want raw data anyway.
# write("commanding raws ...", stderr())
# command_raw = df_sample_sheet %>% 
#     select(sample_id, raw_filename) %>% 
#     rowwise() %>% 
#     mutate(raw_filename_splitted = str_split(raw_filename, ";"),#,
#            source_files = paste0("../output/", sample_id, "/trimmed_reads/", raw_filename_splitted, collapse = " "),
#            target_dir = target_raw) %>% 
#     ungroup() %>% 
#     
#     transmute(command = paste("cp", source_files, target_dir))
# 
# command_raw %>% select(`#!/bin/bash` = command) %>% write_tsv(paste0("output/", batch, "/", batch, "_cp_raw.sh"))
# 


# Copy Consensus data
write("commanding consensus' ...", stderr())
command_consensus = df_sample_sheet %>% 
    select(raw_full_name, consensus_filename) %>% 
    rowwise() %>% 
    mutate(source_file = paste0("../output/", raw_full_name, "/consensus/", raw_full_name, ".fa"),
         #target_dir = target_consensus,
         #target_basename = paste0(consensus_filename),
         target_file = paste0(target_consensus, consensus_filename)) %>% 
    ungroup() %>% 
  
    #transmute(command = paste0("cp ", source_file, " ", target_dir, target_basename))
    transmute(command = paste("cp", source_file, target_file))

# Have a look at the consensus commands
#command_consensus$command[1]

command_consensus %>% select(`#!/bin/bash` = command) %>% write_tsv(paste0("output/", batch, "/", batch, "_cp_consensus.sh")) # TODO: Make as a space-delimited file instead of tsv.

# if (!development_mode) {
#     system(command_consensus %>% pull(command) %>% paste(collapse = "; "))
# }

# Tar the files
# Could also be called from the workflow. But here it is possibly easier.
# TODO: Consider deleting some intermediary files

#write("compressing ...", stderr())
#system(paste0("cd ", target_copy, "; tar -czvf ../../../", file_targz_out, " *"))




if (fail_on_purpose) {
    stop(fail_msg)
    
}

       
# upload the files via FTP

#write("\nDone for now.", stderr())




