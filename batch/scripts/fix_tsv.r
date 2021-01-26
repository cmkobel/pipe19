library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

write("reading input file ...", stderr())
data = read_tsv(args[1]) %>%
    filter(batch != "batch") %>% 
    format_tsv() %>% 
    read_tsv


write(paste(names(data), collapse = " "), stderr())


# If it goes well, the file should be overwritten.
write("overwriting file ...", stderr())
data %>% write_tsv(args[1])


# consider deleting the temporary file.

write("completed fixing ...", stderr())