
devel_mode = F
library(tidyverse, quietly = T)


# rm(list = ls()); devel_mode = T




args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())




integrated_file = args[1]
arg_batch = args[2]
out_p1 = args[3]

if (devel_mode) {
    integrated_file = "~/GenomeDK/clinmicrocore/pipe19/batch/integrated_init.tsv"
    arg_batch = 210126
    out_p1 = "writethefilehere.pdf"
    
}



integrated = read_tsv(integrated_file) %>% 
    filter(!str_detect(raw_sample_name, "^afd")) %>% 
    filter(!str_detect(tolower(raw_sample_name), "^positiv")) %>%
    mutate(totalMissing = if_else(is.na(totalMissing), 29903, totalMissing)) %>% 
    filter(batch == arg_batch)



size =  round( 1.5 + (integrated %>% pull(plate) %>% unique %>% length)*(1))

# p1 
integrated %>% 
    #filter(plate != 9999) %>% 
    group_by(batch, plate) %>% 
    arrange(totalMissing) %>% 
    mutate(rank_totalMissing = row_number(totalMissing)) %>% 
    
    mutate(type = if_else(type == "control/other", "neg. control", type)) %>% 
    
    select(batch, raw_sample_name, type, totalMissing, rank_totalMissing) %>% 
    
    ggplot(aes(rank_totalMissing, totalMissing)) + 
    geom_hline(yintercept = 29903/2, color = "blue", alpha = 0.4) +
    geom_line(color = "grey50") + 
    
    geom_point(aes(color = type)) +
    
    facet_grid(paste0(batch, "\n", plate)~., scales = "fixed") +
    

    xlim(1, 100) +
    ylim(-200, 29903+200) +
    labs(title = paste("Batch", arg_batch),
         x = "samples ordered by missing bases",
         y = "missing bases",
         caption = "The blue line indicates coverage in half a SARS-CoV-2 genome.")


ggsave(out_p1, height = size, width = 10.5)

