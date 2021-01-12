#rm(list = ls())

library(tidyverse, quietly = T)


write("", stderr())

args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

batch = args[1]
plate = args[2]
path = args[3]
year = args[4]
convertRIV = args[5]

# Example call:
# RScript path/to/parse_path.r 210108 20 ~/GenomeDK/ClinicalMicrobio/faststorage/BACKUP/N331/210108_NS500158_0512_AHGLJ7AFX2/fastq TRUE

# For development, the arguments can also be given manually:
#batch = 210108
#plate = 3471
#path = "~/GenomeDK/ClinicalMicrobio/faststorage/BACKUP/N331/210108_NS500158_0512_AHGLJ7AFX2/fastq"
#year = 20
#convertRIV = "TRUE" # You can't pass a type boolean over cli. Only text



files = list.files(path)
write(paste(length(files), "files in", path), stderr())


input = read_table(files, col_names = "basename") %>% 
    mutate(basename_duplicate = basename) %>% 
    separate(basename_duplicate, c("sample_name", "source_project", "moma_serial", "illumina_serial", "lane", "direction", "extension"), "(_|-)") %>% 
    separate(extension, c("001", "extension"), 3) 

# convert 88, 89, 90  ->  R, I, V
if (convertRIV == "TRUE") {
    input = input %>% 
        #separate(sample_name, c("mads_type", "mads_sub_sample_name"), 2) 
        mutate(sample_name_prefix = str_sub(sample_name, 1, 2),
               sample_name_suffix = str_sub(sample_name, 3),
               markRIV = if_else(str_detect(sample_name_prefix, "(88|89|90)"), T, F),
               sample_name_prefix_converted = recode(sample_name_prefix,
                                                     "88" = "R",
                                                     "89" = "I",
                                                     "90" = "V"),
               reconst_sample_name = paste0(sample_name_prefix_converted, sample_name_suffix)) %>% 
        
        # Make it look like it never happened
        select(basename, sample_name = reconst_sample_name, source_project, moma_serial, illumina_serial, lane, direction, `001`, extension)
        
}

 
#insert year
input = input %>% 
    mutate(sample_name = paste0(str_sub(sample_name, 1, 1), year, str_sub(sample_name, 2)))





  #########################
 # Group lanes and pairs #
#########################

# routine function
pastecollapsed = function(x) {
    paste(x, collapse = " ")
}


# Collect the files for each 
input_grouped = input %>% 
    arrange(lane) %>% 
    pivot_wider(id_cols = c(sample_name, source_project, extension, moma_serial),
                names_from = direction,
                values_from = basename,
                values_fn = pastecollapsed) %>% 
    mutate(path = path,
           batch = batch,
           plate = plate) %>% 
    select(batch, plate, moma_serial, sample_name, extension, source_project, path, R1, R2) # Reorder the columns


# Write the input_grouped table to disk so the python-gwf-workflow can be started.

write("writing parsed input list to stdout", stderr())
input_grouped %>% format_tsv() %>% write(stdout())



    
    
