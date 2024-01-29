setwd("~/work/Olink/YKKY0100")
## 加载依赖包
suppressPackageStartupMessages(library(EPARS))
suppressPackageStartupMessages(library(OlinkAnalyze))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggsci))
library(tidyverse)
library(ggrepel)
library(showtext)
library(showtextdb)
showtext_auto()
#### CSF ####
load("~/work/Olink/YKKY0100/Olink_CSF/results_CSF/rendering_cache_Olink/DEG-calculation_721c515ef58d91962b9a95096dd3b72b.RData")
deg <- exp_group_object@deg %>% names() %>% 
  lapply(\(.x) {
    .data <- exp_group_object@deg[[.x]]
    .data@result_df %>% 
     dplyr::filter(Change != "NOT") %>% 
      dplyr::select(Mean_difference, Pvalue) %>% 
      rownames_to_column("Gene") %>% mutate(Group = .x) %>% 
      mutate(label = ifelse(Pvalue < 0.01, "P-val<0.01","P-val >= 0.01"))
  }) %>% Reduce(rbind, .)

# plot backgroud data
bgd <- deg %>% group_by(Group) %>% 
  summarise_all(list(min = min, max = max)) %>% 
  dplyr::select(Group, Mean_difference_min, Mean_difference_max, label_min) %>% 
  rename(label = label_min)
bgd[bgd$Mean_difference_min > 0, ]$Mean_difference_min <- 0
bgd[bgd$Mean_difference_max < 0, ]$Mean_difference_max <- 0
# if min value bigger than 0, the min value is 0
# library(ggrepel)
  
p <- ggplot() +
        geom_col(data = bgd, mapping = aes(x = Group, y = Mean_difference_min),
                 fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
        geom_col(data = bgd, mapping = aes(x = Group, y = Mean_difference_max),
                 fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
        geom_jitter(data = deg, aes(x = Group, y = Mean_difference, color = label),
                    position = position_jitter(seed = 10085), size = 2.5) +
        scale_color_manual(values = c("#db5a6b", "Black")) +
        geom_tile(data = deg, aes(x = Group, y = 0, fill = Group),
                  color = "black", alpha = 0.6, show.legend = F, height = 0.3) +
        geom_text(data = deg, aes(x = Group, y = 0, label = Group), size = 5, 
                  color = "black")+
        ggsci::scale_fill_npg() + # 自定义颜色
        geom_text_repel(data = deg %>% dplyr::filter(abs(Mean_difference) > 0.3), 
                        aes(x = Group, y = Mean_difference, label = Gene),
                        position = position_jitter(seed = 10085), max.overlaps = 13) +
        theme_minimal() + 
        ylab("Mean difference") +
        theme(axis.title = element_text(family = "sans", size = 13,color = "black",face = "bold"),
              axis.text = element_text(family = "sans", colour = "black", size = 12),
          axis.line.y = element_line(color = "black",size = 1),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          # axis.text.y = element_text(family = "sans", colour = "red", size = 9),
          panel.grid = element_blank(),
          text = element_text(family = "sans", colour = "black", size = 9),
          legend.position = "top",
          legend.direction = "vertical",
          legend.justification = c(1,0),
          legend.title = element_blank(),
          legend.text = element_text(size = 13))
ggsave(filename = "./Olink_CSF/multigroup_diff_volcano_plot.pdf", width = 20, height = 12, plot = p)

#### plasma ####
load("~/work/Olink/YKKY0100/Olink_plasma/results_plasma/rendering_cache_Olink/DEG-calculation_55d5c4ba314cbbb7e3ee0dca4a9def16.RData")
deg <- exp_group_object@deg %>% names() %>% 
  lapply(\(.x) {
    .data <- exp_group_object@deg[[.x]]
    .data@result_df %>% 
      dplyr::filter(Change != "NOT") %>% 
      dplyr::select(Mean_difference, Pvalue) %>% 
      rownames_to_column("Gene") %>% mutate(Group = .x) %>% 
      mutate(label = ifelse(Pvalue < 0.01, "P-val<0.01","P-val >= 0.01"))
  }) %>% Reduce(rbind, .)
deg$Group2 <- gsub(x = deg$Group, pattern = "\\(", replacement = " \\(")
# plot backgroud data
bgd <- deg %>% group_by(Group) %>% 
  summarise_all(list(min = min, max = max)) %>% 
  dplyr::select(Group, Mean_difference_min, Mean_difference_max, label_min) %>% 
  rename(label = label_min)
# if min value bigger than 0, the min value is 0
bgd[bgd$Mean_difference_min > 0, ]$Mean_difference_min <- 0
bgd[bgd$Mean_difference_max < 0, ]$Mean_difference_max <- 0
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
p <- ggplot() +
  geom_col(data = bgd, mapping = aes(x = Group, y = Mean_difference_min),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_col(data = bgd, mapping = aes(x = Group, y = Mean_difference_max),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_jitter(data = deg, aes(x = Group, y = Mean_difference, color = label),
              position = position_jitter(seed = 11), size = 2.5) +
  scale_color_manual(values = c("#db5a6b", "Black")) +
  geom_tile(data = deg, aes(x = Group, y = 0, fill = Group),
            color = "black", alpha = 0.6, show.legend = F, height = 0.3) +
  geom_text(data = deg, aes(x = Group, y = 0, label = str_wrap(Group2, width = 1)), size = 4, 
            color = "black") +
  scale_fill_manual(values = col_vector)+
  geom_text_repel(data = deg %>% dplyr::filter(abs(Mean_difference) > 0.3), 
                  aes(x = Group, y = Mean_difference, label = Gene),
                  position = position_jitter(seed = 11), max.overlaps = 21) +
  theme_minimal() + 
  ylab("Mean difference") +
  theme(axis.title = element_text(family = "sans", size = 13,color = "black",face = "bold"),
        axis.text = element_text(family = "sans", colour = "black", size = 12),
        axis.line.y = element_line(color = "black",size = 1),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.y = element_text(family = "sans", colour = "red", size = 9),
        panel.grid = element_blank(),
        text = element_text(family = "sans", colour = "black", size = 9),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.title = element_blank(),
        legend.text = element_text(size = 13))
ggsave(filename = "./Olink_plasma/multigroup_diff_volcano_plot.pdf", width = 20, height = 12, plot = p)



