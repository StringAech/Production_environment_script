immcell_exp <- read.csv("~/work/DSP/YKKY0273/DSP-WTA/results_newgroup/6.Signature/ImmCellAI/ImmCellAI_hsa_result.csv") %>% 
  dplyr::select(-InfiltrationScore)
exp_group_object <- readRDS("~/work/DSP/YKKY0273/DSP-WTA/exp_group_object.RDS")
exp_group_object@deg <- list()
immcell_exp_ratio <- apply(immcell_exp, 1, \(.x) .x/sum(.x) ) %>% t %>% as.data.frame()
anno <- exp_group_object@anno %>% 
  dplyr::select(AOI_Group1)
immcell_exp_ratio <- merge(anno, immcell_exp_ratio, by = 0, sort = F)
immcell_exp_ratio <- immcell_exp_ratio %>% 
  dplyr::select(SampleID = Row.names, everything()) %>% 
  pivot_longer(cols = 3:26, names_to = "CellType", values_to = "Ratio")
## 在X轴添加注释
p <- immcell_exp_ratio %>% 
  ggbarplot(x = 'SampleID', y = 'Ratio', fill = 'CellType') +
  xlab("") +
  scale_fill_manual(values = col_vector) +
  theme(# axis.text.x = element_text(angle = 90,
    # hjust = 1,
    # vjust = 1),
    axis.text.x = element_blank(),
    legend.position = 'right'
  ) + 
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(~AOI_Group1, scales = "free", switch = "x", space = "free") +
  theme(panel.spacing = unit(0, "mm"))