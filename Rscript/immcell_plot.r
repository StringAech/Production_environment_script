library(openxlsx)
library(ggplot2)
library(tidyverse)
library(tidyverse)
suppressPackageStartupMessages(library(RColorBrewer))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


immcell_exp <- read.csv("./results/6.Signature/ImmCellAI/ImmCellAI_hsa_result.csv")
immcell_exp_ratio <- apply(immcell_exp[, 1:24], 1, \(.x) .x/sum(.x) ) %>% t %>% as.data.frame()

immcell_exp_ratio <- immcell_exp_ratio %>% 
  rownames_to_column("SampleID") %>%
  pivot_longer(cols = 2:25, names_to = "CellType", values_to = "Ratio") %>% 
  dplyr::mutate(SampleID = SampleID %>% str_extract("(?<=A\\_).*"), 
                Group = SampleID %>% str_extract("(?<=[0-9]\\_).*")) %>% 
  arrange(Group)

# immcell_exp_ratio$SampleID <- immcell_exp_ratio$SampleID %>% str_extract("(?<=A\\_).*")
p <- immcell_exp_ratio %>% 
  ggboxplot(x = 'CellType', y = 'Ratio', fill = "CellType") +
  xlab("")+
  scale_fill_manual(values = col_vector) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 1),
  legend.position = 'none'
  )   
ggsave(filename = "./results/6.Signature/ImmCellAI/Ratio_boxplot.pdf", plot = p, width = 8, height = 8)
p <- immcell_exp_ratio %>% 
  ggbarplot(x = 'SampleID', y = 'Ratio', fill = 'CellType') +
  xlab("") +
  scale_fill_manual(values = col_vector) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 1),
        legend.position = 'right'
        ) 
ggsave(filename = "./results/6.Signature/ImmCellAI/Ratio_barplot.pdf", plot = p, width = 16, height = 8)
#### diff ####
immcell_exp <- immcell_exp[, 1:24] %>% 
  rownames_to_column("SampleID") %>% 
  dplyr::filter(grepl("Other$", SampleID))
anno <- exp_group_object@anno %>% 
  dplyr::select(starts_with("Group")) %>% 
  rownames_to_column("SampleID") %>% 
  dplyr::filter(grepl("Other$", SampleID))
anno <- list(GroupA = anno[, 1:2],
             GroupB = anno[, c(1, 3)],
             GroupC = anno[, c(1, 4)]
             ) %>% lapply(\(.x) .x %>% dplyr::filter(.[, 2] != "Other"))
source("~/script/Rscript/all.function.R")

lapply(names(anno), function(.x){
  print(.x)
  data <- merge(anno[[.x]], immcell_exp, sort = F)
  write.xlsx(data, file = paste0("./results/immcell_boxplot/", .x, "/immcell_expr.xlsx"))
  data <- cbind(data[, 1:2],
                data[, 3:ncol(data)] %>% .[, apply(., 2, sum)!=0])
  adjusted_boxplot_nojitter(input_df = data, 
                            id_colnames = colnames(data)[1:2], 
                            x_aes = colnames(data)[2], 
                            plot_dir = paste0("./results/immcell_boxplot/", .x, "/"), 
                            plot_suffix = "_Group1-VS-Group2.pdf", 
                            comparisons_list = c("Group1", "Group2"), 
                            x_order = c("Group1", "Group2"),
                            plot_width = 12,
                            ncol = 6)
  # write.xlsx(data, file = paste0("./results/immcell_boxplot/", .x, "/immcell_expr.xlsx"))
})




