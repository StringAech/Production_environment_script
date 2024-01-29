load("~/work/DSP/YKKY0153广东省人民医院梅老师/20231220_change_the_group/DSP_WTA/results_change_logfc/rendering_cache_DSP_WTA/signature_analysis_background_b573548a0ad6e741f8d2dbb1dcd32760.RData")
library(dplyr)
library(limma)
library(tidyverse)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(ggtheme)
library(ggprism)
patway_expr <- exp_group_object@signature_list %>% 
  do.call(rbind, .) %>% rownames_to_column("pathway") %>% 
  separate(pathway, c("Class", "pathway.names"), "\\.")
write.xlsx(patway_expr, file = "./personal_analysis/hallmarke_expr.xlsx")
rm(exp_group_object)
#### calculate the DEG ####

anno <- read.xlsx("personal_analysis/只做ssGSEA和ImmunecellAI分析.xlsx") %>% 
  separate(3, c(LETTERS[1:4]), "_") %>% 
  dplyr::mutate(SampleID = `样本编号（报告上展示）必填`) %>% 
  dplyr::mutate(SampleID = paste0(SampleID, "_", D)) %>% 
  dplyr::select(28, 13:27) %>% 
  dplyr::mutate(`Group1-1` = gsub("-", "_neg", `Group1-1`) %>% gsub("\\+", "_pos", .)) %>% 
  apply(., 2, \(.x) { gsub("\\-", "_", .x)}) %>% 
  apply(., 2, \(.x) { gsub(" ", "_", .x)})
colnames(anno) <- gsub("-", "_", colnames(anno))
anno <- as.data.frame(anno)

celltype <- colnames(data)[2:25]
group <- "Group1-1"


deg_groups <- "Group1_1(EGFR_neg:EGFR_pos); HER2(TLS_N:TLS_P); KRAS_1(TLS_N:TLS_P); EGFR_L858R_1(TLS_N:TLS_P); EGFR_19del_1(TLS_N:TLS_P); ALL_Group2_1(EGFR_N:EGFR_P); ALL_Group3_1(EGFR_N:EGFR_P); EGFR_pos_Group1_1(TLS_N:TLS_P); EGFR_neg_Group1(TLS_N:TLS_P); EGFR.19del对比L858R(EGFR_P_L858R:EGFR_P_19del); TLS_P(B_cell_area:T_cell_area); TLS_N(B_cell_area:T_cell_area); ALL_T_Bcell(B_cell_area:T_cell_area); T_cell(TLS_N:TLS_P); B_cell(TLS_N:TLS_P)"




deg_object_list <- setup_deg_object_from_string(deg_groups)
g <- names(deg_object_list)

Plot.function2 <- function(g, anno, patway_expr){
  
cell.anno <- anno[!is.na(anno[[deg_object_list[[g]]@groupinfo[[1]]]]),
                  c("SampleID",deg_object_list[[g]]@groupinfo[[1]])]

cell.anno[[deg_object_list[[g]]@groupinfo[[1]]]] <- factor(cell.anno[[deg_object_list[[g]]@groupinfo[[1]]]],
                                                           levels = c(deg_object_list[[g]]@groupinfo[[3]],
                                                                      deg_object_list[[g]]@groupinfo[[2]]))
cell.anno <- cell.anno %>% arrange(.[[2]])
lev <- table(cell.anno[, 2])
cal.data <- patway_expr %>% dplyr::select(2, cell.anno[, 1])
rownames(cal.data) <- cal.data[, 1]
cal.data <- cal.data[, 2:ncol(cal.data)]

if(all(colnames(cal.data) == cell.anno$SampleID)){
  
  dir <- paste0("./personal_analysis/hallmarker.boxplot/", g)
  dir.create(dir, recursive = T)
  
  group <- factor(c(rep("Case", as.numeric(lev)[1]),
                    rep("control", as.numeric(lev)[2]) ), 
                  levels = c("Case", "control") )
  
  design <- model.matrix(~0+group)
  colnames(design) <- levels(factor(group))
  rownames(design) <- colnames(cal.data)
  
  compare <- makeContrasts(Case - control, levels=design)
  fit <- lmFit(cal.data, design)
  fit2 <- contrasts.fit(fit, compare)
  fit3 <- eBayes(fit2)
  Diff <- topTable(fit3, coef=1, number=200) #%>% rownames_to_column("PathWay")
  write.xlsx(Diff, file = paste0(dir, "/", g, "_pathwayDiff_limma.xlsx"), overwrite = T, rowNames = T)
  #### plot ####
  dat_plot <- data.frame(id = row.names(Diff),t = Diff$t)
  dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, 
                                     ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),
                              levels=c('Up','Down','NoSignifi'))
  dat_plot <- dat_plot %>% arrange(t)
  # 变成因子类型
  dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
  
  p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
    geom_col()+
    coord_flip() +
    scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
    geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
    xlab('') + 
    ylab('t value of GSVA score, tumour versus non-malignant') + #注意坐标轴旋转了
    guides(fill=F) + # 不显示图例
    theme_prism(border = T) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
      ) + 
    geom_text(data = dplyr::filter(dat_plot, t < -2), 
              aes(x = id,y = 0.1,label = id),
                     hjust = 0,color = 'black') + # 小于-1的为黑色标签
    geom_text(data = dplyr::filter(dat_plot, t < 0 & t > -2),
              aes(x = id,y = 0.1,label = id),
              hjust = 0,color = 'grey') + # 灰色标签
    geom_text(data = dplyr::filter(dat_plot, t > 0 & t < 2),
              aes(x = id,y = -0.1,label = id),
              hjust = 1,color = 'grey') + # 灰色标签
    geom_text(data = dplyr::filter(dat_plot, t >2),
              aes(x = id,y = -0.1,label = id),
              hjust = 1,color = 'black')
  ggsave(plot = p, filename = paste0(dir, "/", g, "_barplot.pdf"), width = 8, height = 8)
  
}else{
  stop("The expression sample and comment sample names do not correspond")
  }

} 



# Plot.function2(g, anno, patway_expr)




for (g in names(deg_object_list)) {
  Plot.function2(g, anno, patway_expr)
}



