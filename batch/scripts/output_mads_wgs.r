library(tidyverse)


args = commandArgs(trailingOnly = TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

# Read cli-args
#batch = args[1]
integrated_file = args[1]
arg_batch = args[2]
file_out = args[3]



devel = F
# devel = T

if (devel) {
    rm(list = ls())
    devel = T
    
    integrated_file = "~/GenomeDK/clinmicrocore/pipe19/batch/integrated.tsv" # Should not be a batch specific file, as it is important to filter for the newest batch.
    
    arg_batch = "210218"
  
    file_out = paste0("~/GenomeDK/clinmicrocore/pipe19/batch/mads/output/32092_WGS_", arg_batch, ".csv")
      
      
}
    




integrated = read_tsv(integrated_file) %>% 
    #select(final_sample_name, raw_full_name, batch, plate, moma_serial, raw_sample_name, ya_sample_name, type, extension, source_project, path, totalMissing, `cprnr.`, koen, Alder, navn, proevenr, afsendt, modtaget, proevekategori, anatomi, ct) %>% 
    
    select(ya_sample_name, batch, type, afsendt, lineage = lineage, clade, totalMissing_interpreted, plate_control_summary) %>% 

    
    
    # Remove potential duplicates
    group_by(ya_sample_name) %>% 
    arrange(desc(batch)) %>% # Put the newest batch first
    mutate(rank = row_number(ya_sample_name)) %>% 
    filter(rank == 1) %>% # Pick the newest. Conclusion: When you later filter for batch, you will have maked sure that each sample name can only be written into mads once.
    ungroup() %>%   
    
    
    # After checking that there are no dupblicates, we can filter by sample
    filter(type == "sample") %>% 
    
    #filter(batch == arg_batch) %>% 
    
    
    # Pangolin marks lineage as "None" when it is missing. I want it to be coded as NA instead
    mutate(lineage = case_when(lineage != "None" ~ lineage))   

    # remove duplicates 



# Throw a warning, if NAs are present in the imported table

nalen = integrated %>% filter(is.na(ya_sample_name)) %>% 
    pull(ya_sample_name) %>% 
    length

if (nalen == 1) {
    write(paste("Warning: imported table contains", nalen, "missing sample names"), stderr())
} else {
    write(paste("Info: all sample names are present in imported data table"), stderr())
}
  



# sseqlist14 = c("I346071", "I367762", "I386326", "I386447", "I386906", "I387313",
#                "I387711", "P000881", "P604057", "R183368", "R185218", "R185725", "R186025", "R186032",
#                "R186033", "R186037", "R186136", "R186141", "R187001", "R187160", "R187422", "R187424",
#                "R187623", "R188503", "R275537", "R275587", "R275763", "R275900", "R275907", "Seqpos", 
#                "V265768")
# 
# sseqlist16 =  c("I386994", "I387783", "I387808", "P604333", "R188714", "R275127",
#                 "R275133", "Seqpos")




# Make the WGS-table
out = integrated %>% 
    filter(!is.na(ya_sample_name)) %>% 
  
  
    filter(batch == arg_batch) %>% 
    

    
    # Call WGS smitsomhed
    mutate(inconclusive = if_else(is.na(lineage) | is.na(clade) | plate_control_summary == "unsatisfactory", T, F), # If lineage or clade is NA, the sample should answer as inconclusive
           MDSU = 32092) %>% 
           
    # Format WGS strings using the inconclusive control for guidance.
    mutate(WGS_smitsomhed = case_when(inconclusive ~ "WGS: Prøven er ikke sekventerbar",
                                      lineage == "B.1.1.7" ~ "WGS: Variant med øget smitsomhed ", # Engelsk variant
                                      lineage == "B.1.351" | lineage == "P.1" ~ "WGS: Variant med øget smitsomhed og nedsat følsomhed for antistoffer",
                                      lineage == "B.1.525" ~ "WGS: Variant med nedsat følsomhed for antistoffer",
                                      TRUE  ~ "Variant med formodet normal smitsomhed"), # else
           WGS_linje = case_when(inconclusive ~ "WGS: inkonklusiv",
                                 TRUE ~paste0("WGS: ", lineage, ", ", clade))) 
    
    # Format for mads-import
    

out %>% 
    select(`sample-id` = ya_sample_name, MDSU, WGS_linje, WGS_smitsomhed) %>% 
    pivot_longer(c(WGS_linje, WGS_smitsomhed)) %>% 
    
    #write_delim(paste0("~/GenomeDK/clinmicrocore/pipe19/batch/mads/output/32092_WGS_", arg_batch, ".csv"), delim = ";")  # TODO: set output path for args
    #arrange(`sample-id`, name) %>% # DO NOT SORT
    
    #write.table(paste0("~/GenomeDK/clinmicrocore/pipe19/batch/mads/output/32092_WGS_", arg_batch, ".csv"), quote = F, sep = ";", fileEncoding = "cp1252", row.names = F)  # TODO: set output path for args
    write.table(file_out, quote = F, sep = ";", fileEncoding = "cp1252", row.names = F)  # TODO: set output path for args







