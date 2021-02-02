
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
out_p2 = args[4]

if (devel_mode) {
    integrated_file = "~/GenomeDK/clinmicrocore/pipe19/batch/integrated_init.tsv"
    arg_batch = 210120
    out_p1 = "writethefilehere_A.pdf"
    out_p2 = "writethefilehere_B.pdf"
    
}



integrated = read_tsv(integrated_file) %>% 
    filter(!str_detect(raw_sample_name, "^afd")) %>% 
    filter(!str_detect(tolower(raw_sample_name), "^positiv")) %>%
    #mutate(totalMissing = if_else((is.na(totalMissing) | totalMissing == "totalMissing"), 29903, as.double(totalMissing))) %>%
    mutate(totalMissing = if_else((is.na(totalMissing) | totalMissing == "totalMissing"), 29903, as.double(totalMissing))) %>%
    
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
    geom_hline(yintercept = 29903/2, color = "blue", alpha = 0.3) +
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


integrated %>% 
    #filter(plate != 9999) %>% 
    group_by(batch, plate) %>% 
    arrange(totalMissing) %>% 
    mutate(rank_totalMissing = row_number(moma_serial)) %>% 
    
    mutate(type = if_else(type == "control/other", "neg. control", type)) %>% 
    
    select(batch, raw_sample_name, type, totalMissing, rank_totalMissing) %>% 
    
    ggplot(aes(rank_totalMissing, totalMissing)) + 
    geom_hline(yintercept = 29903/2, color = "blue", alpha = 0.3) +
    #geom_line(color = "grey50") + 
    
    geom_point(aes(color = type)) +
    
    facet_grid(paste0(plate, "\n", plate_negative_control_summary)~., scales = "fixed") +

    
    scale_x_continuous(breaks = seq(0, 96, 8), minor_breaks = NULL)+
    scale_y_continuous(breaks = c(0, 29903), limits = c(-10000, 40000))+ 


    labs(title = paste("Batch", arg_batch),
         x = "samples ordered by library ID (moma_serial)",
         y = "missing bases",
         caption = "The blue line indicates coverage in half a SARS-CoV-2 genome.")



ggsave(out_p2, height = size, width = 10.5)


