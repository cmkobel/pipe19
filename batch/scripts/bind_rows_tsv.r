library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())


# random_tmp = paste0("tmp_", sample(1000000:9999999, 1), ".tsv")
# command = paste("cat", paste(args, collapse = " "), ">", random_tmp)
# write(command, stderr())
# system(command)


col_types_string = paste(rep("c", args[1]))
write(paste(col_types_string, collapse = ", "), stderr())

data = tibble()
for (input_table in args[-1]) {
    data = bind_rows(read_tsv(input_table), col_types = col_types_string)
}

#data = read_tsv(random_tmp) %>% filter(sample != "sample") %>% format_tsv



write(data %>% format_tsv, stdout())



# consider deleting the temporary file.