devel = F
#rm(list = ls()); devel = T





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
if (devel) {
    
    batch = 210108
    #plate = 3471
    path = "~/GenomeDK/ClinicalMicrobio/faststorage/BACKUP/N331/210108_NS500158_0512_AHGLJ7AFX2/fastq"
    year = 20
    convertRIV = "TRUE" # You can't pass a type boolean over cli. Only text
    format_specifier = "formatA"
}


if (format_specifier == "formatA") {
    column_format = c("sample_name", "source_project", "moma_serial", "illumina_serial", "lane", "direction", "extension")
} else if (format_specifier == "formatB") {
    column_format = c("sample_name", "plate", "source_project", "moma_serial", "illumina_serial", "lane", "direction", "extension")
}


write("reading list of files ...", stderr())
files = list.files(path)
write(paste(length(files), "files in", path), stderr())

if (length(files) < 1) {
  cat("\n\n")
  stop(paste("the path specified does not contain any files. \nPlease check that the path is correct."))
}

write("parsing list of files ...", stderr())
input = read_table(paste0(files, collapse = "\n"), col_names = "basename")%>%  
    mutate(basename_duplicate = basename) %>% 
    separate(basename_duplicate, column_format, "(_|-)") %>% 
    separate(extension, c("001", "extension"), 3) %>% 
    
    # add a column that tells whether the file is a control or not
    mutate(type = if_else(str_detect(tolower(sample_name), "negativ|positiv"),
                          "control",
                          "sample"))


# Backwards compatibility:
# If no plate number is given in the file name, insert 0000
if (format_specifier == "formatA") {
    input = input %>% 
        mutate(plate = "0000")
} # else: do nothing - the plate number will be given by the filename.

 
# Insert the year for the samples
# TODO: This should be a list instead
year_R = paste0("R", year)
year_I = paste0("I", year)
year_V = paste0("V", year)

# convert 88, 89, 90  ->  R, I, V
if (convertRIV == "TRUE") {
write("converting RIV ...", stderr())
    input = input %>% 
        #separate(sample_name, c("mads_type", "mads_sub_sample_name"), 2) 
        mutate(sample_name_prefix = str_sub(sample_name, 1, 2),
               sample_name_suffix = str_sub(sample_name, 3),
               markRIV = if_else(str_detect(sample_name_prefix, "(88|89|90)"), T, F),
               sample_name_prefix_converted = recode(sample_name_prefix,
                                                     "88" = year_R,
                                                     "89" = year_I,
                                                     "90" = year_V),
               reconst_sample_name = paste0(sample_name_prefix_converted, sample_name_suffix)) %>% 
        
        # Make it look like it never happened
        select(basename, plate, sample_name = reconst_sample_name, source_project, moma_serial, illumina_serial, lane, direction, `001`, extension, type)
        
}

 

sn_unique_lengths = input %>% filter(type == "sample") %>% 
    mutate(sn_len = str_length(sample_name)) %>% 
    pull(sn_len) %>% 
    unique


if (length(sn_unique_lengths) != 1 | sn_unique_lengths != 9) {
    
    sn_unique_lengths_pasted = paste(sn_unique_lengths, collapse = " ")
    # TODO: Consider that this should maybe be a stop instead of print
    write(paste("Warning: The lengths of the sample names is", sn_unique_lengths_pasted, ". It ought to be 9 <LYYNNNNNN>"), stderr())
    

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
    select(batch, plate, moma_serial, sample_name, type, extension, source_project, path, R1, R2) # Reorder the columns


# Write the input_grouped table to disk so the python-gwf-workflow can be started.

write("writing parsed input list to stdout", stderr())
input_grouped %>% format_tsv() %>% write(stdout())



    
    
