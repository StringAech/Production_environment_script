setwd("~/work/DSP/YKKY0153广东省人民医院梅老师/20231220_change_the_group/DSP_WTA")
library(openxlsx)
library(dplyr)
library(limma)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)


anno <- read.xlsx("personal_analysis/只做ssGSEA和ImmunecellAI分析.xlsx") %>% 
  separate(3, c(LETTERS[1:4]), "_") %>% 
  dplyr::mutate(SampleID = `样本编号（报告上展示）必填`) %>% 
  dplyr::mutate(SampleID = paste0(SampleID, "_", D)) %>% 
  dplyr::select(28, 13:27)


data <- read.csv("./results/5.Signature/ImmCellAI/ImmCellAI_hsa_result.csv") %>% 
  rownames_to_column("SampleID")
celltype <- colnames(data)[2:25]

group <- "Group1-1"
cell.anno <- anno[!is.na(anno[[group]]), c("SampleID", group)]
cell.data <- merge(data, cell.anno, sort = T) %>% dplyr::select(1, 27, all_of(celltype)) %>% 
  dplyr::select(1, Group = 2, everything())


boxplt.immcell <- function(cell.data, celltype = celltype, dir = "./", title){
 plot.data <- pivot_longer(cell.data, cols = -colnames(cell.data)[1:2], 
                          names_to = "Celltype", values_to = "Abundance")
plot.data$Celltype <- factor(plot.data$Celltype, levels = celltype)
stat.res <- plot.data %>%
  group_by(Celltype) %>%
  rstatix::t_test(Abundance ~ Group, detailed = T) %>%
  add_xy_position(x = "Celltype", scales = "free") %>%
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance("p") %>% 
  dplyr::mutate(y.position = max(plot.data$Abundance) + 0.1)
write.xlsx(stat.res, file = paste0(dir, "/stat.res.xlsx"), overwrite = T)

p <- ggboxplot(plot.data, x = "Celltype", y = "Abundance", fill = "Group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "right") +
  stat_pvalue_manual(stat.res, 
                     label = "p.signif",
                     # label = "p", 
                     tip.length = 0) + 
  ggsci::scale_fill_nejm() + 
  ggtitle(title)
ggsave(filename = paste0(dir, "/immcellabundance_boxplot.pdf"), plot = p, 
      height = 9, width = 16
       )
}

lapply(colnames(anno)[-1], function(.x){
  dir <- paste0("./personal_analysis/immcellAI_change_p/", .x)
  dir.create(dir, recursive = T)
  cell.anno <- anno[!is.na(anno[[.x]]), c("SampleID", .x)]
  cell.data <- merge(data, cell.anno, sort = T) %>% dplyr::select(1, 27, all_of(celltype)) %>% 
    dplyr::select(1, Group = 2, everything())
  boxplt.immcell(cell.data = cell.data, 
                 celltype = celltype, 
                 dir = dir, 
                 title = paste0(.x, "_immcellabundance"))
})






