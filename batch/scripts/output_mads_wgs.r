


library(tidyverse)

rm(list = ls())

arg_batch = "test" # should be set from args



integrated = read_tsv("~/GenomeDK/clinmicrocore/pipe19/batch/integrated.tsv") %>% 
    #select(final_sample_name, raw_full_name, batch, plate, moma_serial, raw_sample_name, ya_sample_name, type, extension, source_project, path, totalMissing, `cprnr.`, koen, Alder, navn, proevenr, afsendt, modtaget, proevekategori, anatomi, ct) %>% 
    
    select(ya_sample_name, batch, type, afsendt, lineage = lineage, clade, totalMissing_interpreted, plate_control_summary) %>% 

    
    
    # Remove potential duplicates
    group_by(ya_sample_name) %>% 
    arrange(desc(batch)) %>% # Put the newest batch first
    mutate(rank = row_number(ya_sample_name)) %>% 
    filter(rank == 1) %>% # Pick the newest 
    ungroup() %>%     
    
    
    # After checking that there are no dupblicates, we can filter by sample
    filter(type == "sample") %>% 
    
    #filter(batch == arg_batch) %>% 
    
    
    # Pangolin marks lineage as "None" when it is missing. I want it to be coded as NA instead
    mutate(lineage = case_when(lineage != "None" ~ lineage))   

    # remove duplicates 



# Throw a warning, if NAs are present in the imported table
{
    nalen = integrated %>% filter(is.na(ya_sample_name)) %>% 
        pull(ya_sample_name) %>% 
        length
    
    if (nalen == 1) {
        write(paste("Warning: imported table contains", nalen, "missing sample names"), stderr())
    } else {
        write(paste("Info: all sample names are present in imported data table"), stderr())
    }
    }



sseqlist14 = c("I346071", "I367762", "I386326", "I386447", "I386906", "I387313",
               "I387711", "P000881", "P604057", "R183368", "R185218", "R185725", "R186025", "R186032",
               "R186033", "R186037", "R186136", "R186141", "R187001", "R187160", "R187422", "R187424",
               "R187623", "R188503", "R275537", "R275587", "R275763", "R275900", "R275907", "Seqpos", 
               "V265768")


# TODO: remove this block for production
# Let's take a few samples only
integrated = integrated %>% 
    arrange(desc(afsendt)) %>% 
    mutate(rn = row_number(afsendt)) %>% 
    #filter(rn <= 5 | ya_sample_name == "I346261" | ya_sample_name == "R138768" | ya_sample_name == "I376722" | ya_sample_name == "I387808") %>% 
    filter(ya_sample_name %in% sseqlist14) %>% 
    select(-rn)


# Make the WGS-table
out = integrated %>% 
    filter(!is.na(ya_sample_name)) %>% 
    

    
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
    arrange(`sample-id`, name) %>% 
    write.table(paste0("~/GenomeDK/clinmicrocore/pipe19/batch/mads/output/32092_WGS_", arg_batch, ".csv"), sep = ";", fileEncoding = "cp1252", row.names = F)  # TODO: set output path for args







