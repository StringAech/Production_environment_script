pathways_difference <- function(pathway_expr, anno, deg_object_list, deg_names=1, dir){
  ## pathway_expr rownames is smapleID, colnames is pathwaynames
  ## anno <- exp_group_object@anno %>% rownames_to_column("SampleID")
  ## deg_object_list is DEG's cache deg_object_list
  ## deg_names is names(exp_group_object@deg)
  ## dir is dir
  # browser()
  groupinfo_col <- deg_object_list[[deg_names]]@groupinfo[[1]]
  cell.anno <- anno %>% 
    dplyr::select(!!sym(groupinfo_col), SampleID) %>% 
    dplyr::filter(!!sym(groupinfo_col) != "Other")
  cell.anno[[groupinfo_col]] <- factor(cell.anno[[groupinfo_col]],
                                       levels = c(deg_object_list[[deg_names]]@groupinfo[[3]],
                                                  deg_object_list[[deg_names]]@groupinfo[[2]]))
  cell.anno <- cell.anno %>% arrange(.[[groupinfo_col]])
  lev <- table(cell.anno[, groupinfo_col])
  cal.data <- pathway_expr %>% 
    as.data.frame() %>% 
    dplyr::select(cell.anno[, "SampleID"])
  if (all(colnames(cal.data) == cell.anno$SampleID)) {
    dir.create(dir, recursive = T, showWarnings = F)
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
    Diff <- topTable(fit3, coef=1, number=200)
    file.name <- paste0(groupinfo_col, 
                        "(", deg_object_list[[deg_names]]@groupinfo[[3]], 
                        "-VS-", deg_object_list[[deg_names]]@groupinfo[[2]], ")")
    write.xlsx(Diff, file = paste0(dir, "/", file.name, "_pathwayDiff_limma.xlsx"), overwrite = T, rowNames = T)
    
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
      geom_hline(yintercept = c(-2, 2), color = 'white', linewidth = 0.5, lty='dashed') +
      xlab('') + 
      ylab(paste0('t value of GSVA score, ', names(lev)[1], " versus ", names(lev)[2])) +
      # ylab('t value of GSVA score, tumour versus non-malignant') + #注意坐标轴旋转了
      guides(fill = "none") + # 不显示图例
      theme_prism(border = T) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) + 
      geom_text(data = dplyr::filter(dat_plot, t < -2), size = 3,# fix the font size
                aes(x = id,y = 0.1, label = id),
                hjust = 0,color = 'black') + # 小于-1的为黑色标签
      geom_text(data = dplyr::filter(dat_plot, t < 0 & t > -2), size = 3,
                aes(x = id,y = 0.1,label = id),
                hjust = 0,color = 'grey') + # 灰色标签
      geom_text(data = dplyr::filter(dat_plot, t > 0 & t < 2), size = 3,
                aes(x = id,y = -0.1,label = id),
                hjust = 1,color = 'grey') + # 灰色标签
      geom_text(data = dplyr::filter(dat_plot, t >2), size = 3,
                aes(x = id,y = -0.1,label = id),
                hjust = 1,color = 'black')
    ggsave(plot = p, filename = paste0(dir, "/", file.name, "_barplot.pdf"), width = 8, height = 12)
  }else{
    stop("The expression sample and comment sample names do not correspond")
  }
}