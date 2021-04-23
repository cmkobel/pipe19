





library(tidyverse, quietly = T)


write("", stderr())

args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

batch = args[1]
#plate = args[2]
path = args[2]
year = args[3]
convertRIV = args[4]
format_specifier = args[5]

# Example call:
# RScript path/to/parse_path.r 210108 20 ~/GenomeDK/ClinicalMicrobio/faststorage/BACKUP/N331/210108_NS500158_0512_AHGLJ7AFX2/fastq TRUE

# For development, the arguments can also be given manually:
if (F) {

    rm(list = ls())
    setwd("~/GenomeDK/clinmicrocore/pipe19/")    
    
    #batch = 210212
    batch = 210421
    #plate = 3471
    #path = "~/GenomeDK/clinmicrocore/BACKUP/N350/210212_NS500158_0522_AHLHF3AFX2/fastq/"
    path = "~/GenomeDK/clinmicrocore/BACKUP/N331/210421_NB501001_0534_AHV27LAFX2/fastq"
    year = 21
    convertRIV = "FALSE" # You can't pass a type boolean over cli. Only text
    format_specifier = "formatB"
}

# Segregate whether the plate id can be found in the filen 
if (format_specifier == "formatA" | format_specifier == "formatC") { 
  
    column_format = c("sample_name", "source_project", "moma_serial", "illumina_serial", "lane", "direction", "extension")
    
    
} else if (format_specifier == "formatC") { # | format_specifier == "formatD") { <-- af en eller anden grund har der også været en D engang.
  
    column_format = c("sample_name", "plate", "source_project", "moma_serial", "illumina_serial", "lane", "direction", "extension")
    
    
# Efter kontraktbruddet var jeg nødt til at skubbe den gamle B til C (ovenfor) og lave dette nye B-format som indeholder "moma_junk"
} else if (format_specifier == "formatB") { 
  
    column_format = c("sample_name", "plate", "moma_junk", "source_project", "moma_serial", "illumina_serial", "lane", "direction", "extension")
  
} else {
  
  stop(paste("format specifier", column_format, "is not supported."))
  
}
# Remember to also update the if statement about formats further down.


write("reading list of files ...", stderr())
files = list.files(path)
write(paste(length(files), "files in", path), stderr())

if (length(files) < 1) {
  cat("\n\n")
  stop(paste("the path specified does not contain any files. \nPlease check that the path is correct."))
}

write("parsing list of files ...", stderr())
input = read_table(paste0(files, collapse = "\n"), col_names = "basename") %>% # write_tsv("debug0.tsv")
    #stop()
    mutate(basename_duplicate = basename) %>% 
    separate(basename_duplicate, column_format, "(_|-)", convert = F) %>% 
    separate(extension, c("001", "extension"), 3, convert = F) %>% #View#write_tsv("debug0.tsv")
    
    # add a column that tells whether the file is a control or not
    # mutate(type = if_else(str_detect(tolower(sample_name), "negativ|positiv|blank|tom|^afd|^00|^h2o|^neg|neg$|^empty"),
    #                       "control/other",
    #                       "sample"))


    mutate(type = case_when(str_detect(tolower(sample_name), "positiv|pos$|seqpos") ~ "positive_control",
                            str_detect(tolower(sample_name), "negativ|h2o|^empty|blank|tom|^neg|neg$") ~ "negative_control",
                            str_detect(tolower(sample_name), "^afd|^00") ~ "other",
                            TRUE ~ "sample"))


write(paste("got to the other side ", getwd()), stderr())
input %>% write_tsv("debug1.tsv")

# Backwards compatibility:
# If no plate number is given in the file name, insert 9999
if (format_specifier == "formatA" | format_specifier == "formatC") {
    input = input %>% 
        mutate(plate = "9999")
} # else: do nothing - the plate number will already have been given by the filename itself.




sn_unique_lengths = input %>% filter(type == "sample") %>% 
    mutate(sn_len = str_length(sample_name)) %>% 
    pull(sn_len) %>% 
    unique

if (FALSE) {
    if (length(sn_unique_lengths) != 1 | sn_unique_lengths != 9) {
        
        sn_unique_lengths_pasted = paste(sn_unique_lengths, collapse = " ")
        # TODO: Consider that this should maybe be a stop instead of print
        write(paste("Warning: The lengths of the sample names is", sn_unique_lengths_pasted, ". It ought to be 9 <LYYNNNNNN>"), stderr())
        
    
    }
}

  #########################
 # Group lanes and pairs #
#########################

# routine function
pastecollapsed = function(x) {
    paste(x, collapse = " ")
}


# Collect the files for each 
write("collecting files for each sample", stderr())
input_grouped = input %>% 
    arrange(lane) %>% 
    pivot_wider(id_cols = c(sample_name, plate, source_project, extension, moma_serial, type),
                names_from = direction,
                values_from = basename,
                values_fn = pastecollapsed) %>% 
    mutate(path = path,
           batch = batch) %>% 
    select(batch, plate, moma_serial, raw_sample_name = sample_name, type, extension, source_project, path, R1, R2) %>% # Reorder the columns
    #arrange(moma_serial)
    arrange(moma_serial)
    


# Write the input_grouped table to disk so the python-gwf-workflow can be started.

write("writing parsed input list to stdout", stderr())
input_grouped %>% format_tsv() %>% write(stdout())



    
    
