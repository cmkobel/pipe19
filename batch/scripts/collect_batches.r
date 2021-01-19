library(tidyverse)

# rm(list = ls()); setwd("~/GenomeDK/ClinicalMicrobio/faststorage/pipe19/batch"); development_mode = T

args = commandArgs(trailingOnly = TRUE)
write(paste("arguments given: ", paste(args, collapse = ", ")), stderr())


out_file = args[1]

data = tibble()
for (batch_integrated in Sys.glob("output/*/*.tsv")) {
    tmp = read_tsv(batch_integrated)
    data = bind_rows(data, tmp)
}



data %>% write_tsv(paste0("backup_all_batches_integrated/", format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), ".tsv"))
data %>% write_tsv(out_fsa)