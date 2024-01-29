exp_group_object@anno %>% colnames()
exp_group_object@deg_groups
library(openxlsx)
read_CTA_anno <- function(anno_path){
  anno_df <- read.xlsx(file.path(anno_path), sheet = 1) %>% as.data.frame()
  ncol_modify_list <- intersect(colnames(anno_df), c("AOISurfaceArea", "AOINucleiCount", "ROICoordinateX", "ROICoordinateY",
                                                     "RawReads", "AlignedReads", "DeduplicatedReads", "TrimmedReads",
                                                     "StitchedReads", "SequencingSaturation", "SequencingSetID", "UMIQ30",
                                                     "RTSQ30", "LOQ.(Cancer.Transcriptome.Atlas)", "ScalingFactor", "NormalizationFactor"))
  anno_df[,ncol_modify_list]<-lapply(anno_df[,ncol_modify_list],as.numeric)
  rownames(anno_df) <- paste(anno_df$ScanLabel,anno_df$ROILabel, anno_df$SegmentLabel, sep="_")
  anno_df$SegmentDisplayName = gsub(" | ", "_", anno_df$SegmentDisplayName, fixed = T, perl = F)
  anno_df$SegmentDisplayName = gsub(" ", ".", anno_df$SegmentDisplayName, fixed = T, perl = F)
  rownames(anno_df) <- gsub(" ", ".", rownames(anno_df))
  anno_df <- anno_df[order(rownames(anno_df)),]
  # anno_mtx <- as.matrix(anno_df)
  return(anno_df)
}

# exp_group_object@anno <- read.xlsx("../YKKY0112_BPQC.xlsx", sheet = 1)

exp_group_object@anno <- read_CTA_anno("../YKKY0112_BPQC.xlsx")

exp_group_object@deg_groups <- c("Morphology(VA:PA),Group$Non_ICI; Morphology(VA:LPA),Group$Non_ICI; Morphology(PA:LPA),Group$Non_ICI; Morphology(VA:PA),Group$ICI; Morphology(VA:LPA),Group$ICI; Morphology(PA:LPA),Group$ICI; Group(Non_ICI:ICI),Morphology2$VAPA; Group(Non_ICI:ICI),Morphology$LPA; Group(Non_ICI:ICI),Morphology$PA; Group(Non_ICI:ICI),Morphology$VA")
CachePath <- "./results/Personalized_analysis/20240119/DEG_data.frame/"
deg_object_list <- setup_deg_object_from_string(exp_group_object@deg_groups)
exp_group_object@deg <- list()
library(parallel);#加载并行计算包
cl <- makeCluster(10);# 初始化cpu集群
do.data=function(g_complex=1){
  # 使用整理好的过滤信息对AOI进行过滤
  count = exp_group_object@exp
  groupinfo <- factor(data.frame(exp_group_object@anno)[[deg_object_list[[g_complex]]@groupinfo[[1]] ]])
  if (deg_object_list[[g_complex]]@filtering_boolean == FALSE && length(deg_object_list[[g_complex]]@groupinfo) == 1){
    deg_object_list[[g_complex]]@filtered_sample <- colnames(exp_group_object@exp)
    deg_object_list[[g_complex]]@filtered_matrix <- count
  } else{
    tmp_filtered_list <- list()
    # 针对每一个AOI筛选条件对每个AOI进行判断
    for (per_filter in 1:length(deg_object_list[[g_complex]]@filtering_list)){
      tmp_filtered_list[[per_filter]] <- eval(parse(text = 
                                                      deg_object_list[[g_complex]]@filtering_list[[per_filter]]))
    }
    # 取交集，所有筛选条件均通过的AOI/ROI被留存下来
    tmp_filtered_list <- Reduce("&",tmp_filtered_list)
    deg_object_list[[g_complex]]@filtered_sample <- colnames(exp_group_object@exp)[tmp_filtered_list]
    deg_object_list[[g_complex]]@filtered_matrix <- count[ ,tmp_filtered_list]
    groupinfo <- factor(groupinfo[tmp_filtered_list])
  }
  # system(paste("echo",g_complex,">>run.txt"))
  
  #输入差异分组时,要求control在前,case在后
  .group <- deg_object_list[[g_complex]]@groupinfo[[1]]
  .case <- deg_object_list[[g_complex]]@groupinfo[[3]]
  .control <- deg_object_list[[g_complex]]@groupinfo[[2]]
  #按照输入的分组顺序进行排序
  if(length(unique(groupinfo))>1){
    groupinfo <- forcats::fct_relevel(groupinfo,c(.control,.case))
    # prefix
    prefix_string <- stringr::str_extract(deg_object_list[[g_complex]]@label_string,".*中")
    label_string <- paste0(prefix_string,.case,"比",.control,"的差异")
    deg_object_list[[g_complex]]@label_string <- label_string
    tmp_deg <- obtain_edgeR_DE(matrix = deg_object_list[[g_complex]]@filtered_matrix,
                               groupinfo = groupinfo, 
                               lfc_cutoff = exp_group_object@config$deg$lfc_cutoff,
                               fdr_cutoff = 1, 
                               outdir = CachePath, 
                               prefix = paste0("DE_list_" ,deg_object_list[[g_complex]]@label_string))
    system(paste("echo",g_complex,">>run.txt"))
    tmp_deg=data.frame(tmp_deg,check.names = F)
    data=new('deg_data')
    #data@groupinfo
    data@result_df=tmp_deg
    data@filtered_matrix=deg_object_list[[g_complex]]@filtered_matrix
    data@filtered_sample=deg_object_list[[g_complex]]@filtered_matrix%>%colnames()
    data@groupinfo=deg_object_list[[g_complex]]@groupinfo
    prefix_string <- stringr::str_extract(deg_object_list[[g_complex]]@label_string,".*中")
    label_string <- paste0(prefix_string,.case,"比",.control,"的差异")
    data@label_string <- label_string
    
    # save the deg result back into the deg_object_list, in the result_df
    return(data)
  }
  system(paste("echo",g_complex,">>run.txt"))
  return(data.frame())
}
clusterExport(cl,'exp_group_object');
clusterExport(cl,'deg_object_list');
clusterExport(cl,'CachePath');
clusterEvalQ(cl,library(dplyr));#添加并行计算中用到的包
clusterEvalQ(cl,library(EPARS));#添加并行计算中用到的包
get1=parLapply(cl,1:length(deg_object_list), do.data)
for (g_complex in 1:length(deg_object_list)){
  deg_object_list[[g_complex]] <- get1[[g_complex]]    
}
parallel::stopCluster(cl)
for (i in names(deg_object_list)){
  exp_group_object@deg[[i]] <- deg_object_list[[i]]
}
for (g in names(exp_group_object@deg)){
  tmp_groupinfo <- exp_group_object@deg[[g]]@groupinfo[[1]]
  if (length(exp_group_object@deg[[g]]@groupinfo) >= 2){
    tmp_groupinfo_1 <- exp_group_object@deg[[g]]@groupinfo[[2]]
    tmp_groupinfo_2 <- exp_group_object@deg[[g]]@groupinfo[[3]]
  } else{
    tmp_groupinfo_1 <- levels(factor(exp_group_object@anno[[exp_group_object@deg[[g]]@groupinfo[[1]]]]))[1]
    tmp_groupinfo_2 <- levels(factor(exp_group_object@anno[[exp_group_object@deg[[g]]@groupinfo[[1]]]]))[2]
  }
  
  exp_group_object@deg[[g]]@natural_annotation <- paste0(sprintf("\n\n 此差异分析研究的是从 **%s** 的数值中分析 **%s** 和 **%s** 两个组别的AOI/ROI 之间的差异表达。如果一个基因在此差异表达中**上调(UP)**，代表着 **%s** 的AOI/ROI中此基因的表达相对 **%s** 的AOI/ROI表达**增高**，相反**下调(DOWN)**代表 **%s** 的AOI/ROI中此基因的表达相对 **%s** 的AOI/ROI表达**降低**。",   paste0(tmp_groupinfo),
                                                                 paste0(tmp_groupinfo_1), paste0(tmp_groupinfo_2),
                                                                 paste0(tmp_groupinfo_2), paste0(tmp_groupinfo_1),
                                                                 paste0(tmp_groupinfo_2), paste0(tmp_groupinfo_1)))
}
saveRDS(exp_group_object, file = "./results/Personalized_analysis/20240119/exp_group_object.RDs")
######
PlotsPath <- "./results/Personalized_analysis/20240119/spatial-circ-boxplot/"
num_AOI <- nrow(exp_group_object@anno)
if("SlideName" %in% names(exp_group_object@anno)){
  unique_slide <- unique(exp_group_object@anno$SlideName)
  num_slide <- length(unique_slide)
}

# 循环title
template <- "#### 差异分组: %s {.unlisted .unnumbered}
"
title_template <- "### 玻片: %s {.tabset .tabset-fade .unnumbered}
"

object <- exp_group_object

# 外圈循环， 以玻片为单位
for (per_slide in unique_slide){
  cat(sprintf(title_template, paste0(per_slide)))
  
  # 内圈循环，以差异分析为单位
  for (g in names(object@deg)){
    cat(sprintf(template, paste0(object@deg[[g]]@label_string)))
    
    # 提取差异分析中最显著的基因个数
    sig_gene_list <- rownames(subset(object@deg[[g]]@result_df, 
                                     change != 'NOT') %>% slice_min(PValue, n = 10))
    # 提取本玻片上的AOI队列
    AOI_list <- rownames(subset(object@anno, SlideName == per_slide))
    
    # 筛选差异分析分组包含的AOI
    AOI_list<- AOI_list[AOI_list %in% object@deg[[g]]@filtered_sample]
    
    # 开始绘图模块
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(plotrix))
    suppressPackageStartupMessages(library(ggforce))
    
    # 提取坐标
    coordinates <- t(object@anno[,c("ROICoordinateX", "ROICoordinateY")])
    # colnames(coordinates) <- colnames(coordinates) %>% str_extract("(?<=WTA\\_).*")
    # 筛选矩阵只保留显著10个基因
    exp <- object@exp[rownames(object@exp) %in% sig_gene_list, ]
    # colnames(exp) <- colnames(exp) %>% str_extract("(?<=WTA\\_).*")
    plotting_data <- as.data.frame(t(rbind(coordinates, exp)))
    # 过滤样本只保留本玻片的样本
    plotting_data <- plotting_data[rownames(plotting_data) %in% AOI_list, ]
    plotting_data$AOI <- rownames(plotting_data)
    plotting_data <- plotting_data %>% tidyr::pivot_longer(-c(ROICoordinateX,ROICoordinateY, AOI))
    plotting_data <- plotting_data %>% group_by(AOI) %>% 
      mutate(start_location = (1:n()-1) /n()*pi*2, 
             end_location = 1:n() /n()*pi*2)
    check_overplot <- as.matrix(unique(data.frame(
      ROICoordinateX = plotting_data$ROICoordinateX,
      ROICoordinateY = plotting_data$ROICoordinateY)))
    class(check_overplot) <- "numeric"
    plotting_data <- data.frame(check_overplot) %>% full_join(plotting_data)
    plotting_data$AOI <- plotting_data$AOI %>% str_extract("(?<=WTA\\_).*")
    p <- ggplot(plotting_data) + geom_arc_bar(aes(x0 = ROICoordinateX, y0 = -ROICoordinateY,
                                                  r0 = 0, r = sqrt(value * 8e7 / max(value)), start = start_location,
                                                  end   = end_location, fill  = name), alpha = 0.65, color = "transparent") + 
      # geom_text(aes(x = ROICoordinateX, y = -ROICoordinateY, label = AOI), size = 2,
      #           hjust = 0.5, vjust = 3, check_overlap = T) +
      coord_fixed() + theme_bw() +
      scale_fill_manual("Genes", values = object@config$pan_color) +
      xlab("Slide_Coordinate_X") + ylab("Slide_Coordinate_Y") +
      theme(panel.grid=element_blank())
    subchunkify(p, fig_height = 15, fig_width = 12)
    ggsave(paste0(PlotsPath, "Spatial_expression_circular_barplot_", gsub(" ", "_", g),
                  "_on_slide_", per_slide, ".pdf"), plot = p, height = 15, width = 12)
  }
}


#### 3 ####
library(stringi)
library(tidyverse)
genelist <- read.xlsx("../Gene list.xlsx")
VAsample <- exp_group_object@anno %>% 
  dplyr::filter(Morphology == "VA") %>% rownames()
PAsample <- exp_group_object@anno %>% 
  dplyr::filter(Morphology == "PA") %>% rownames()
PLAsample <- exp_group_object@anno %>% 
  dplyr::filter(Morphology == "LPA") %>% rownames()
colnames(exp_group_object@exp) <- colnames(exp_group_object@exp) %>% 
  str_extract("(?<=\\)\\-).*")
exp_group_object@anno$SegmentDisplayName <- exp_group_object@anno$SegmentDisplayName %>% 
  str_extract("(?<=\\)\\-).*")

Target_genes_expr <- exp_group_object@exp %>% as.data.frame() %>% 
  dplyr::filter(rownames(.) %in% genelist$Target.genes) %>% 
  t %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  mutate(Group = stri_replace_all_regex(SampleID, 
                                        pattern = exp_group_object@anno$SegmentDisplayName,
                                        replacement = exp_group_object@anno$Morphology, vectorize_all = F)) %>% 
  dplyr::select(SampleID, Group, everything())
source("~/script/Rscript/all.function.R")
Target_genes_expr %>%
  adjusted_boxplot_nojitter(input_df = .,
                            id_colnames = colnames(.)[1:2],
                            x_aes = colnames(.)[2],
                            plot_dir = "./results/Personalized_analysis/20240119/boxplot/",
                            plot_suffix = "-VA-VS-LPA.pdf",
                            comparisons_list = list(c("VA", "LPA")),
                            x_order = c("VA", "LPA"),plot_width = 25,ncol = 6
                            )

plotdata.function <- function(genelist, data = exp_group_object@exp, anno = exp_group_object@anno){
  data %<>% as.data.frame() %>% 
    dplyr::filter(rownames(.) %in% genelist) %>% 
    t %>% 
    as.data.frame() %>% 
    rownames_to_column("SampleID") %>% 
    mutate(Group = stri_replace_all_regex(SampleID, 
                                          pattern = anno$SegmentDisplayName,
                                          replacement = anno$Morphology, vectorize_all = F)) %>% 
    dplyr::select(SampleID, Group, everything())
  return(data)
}
listgene <- list()
for (i in colnames(genelist)){
  listgene[[i]] <- genelist[, i]
}
expr_list <- lapply(listgene, function(.x){
  .data <- plotdata.function(genelist = .x, data = exp_group_object@exp, anno = exp_group_object@anno)
  return(.data)
})

lapply(names(expr_list), function(.x){
  item <- expr_list[[.x]][["Group"]] %>% unique()
  expr_list[[.x]] %>%
    adjusted_boxplot_nojitter(input_df = .,
                              id_colnames = colnames(.)[1:2],
                              x_aes = colnames(.)[2],
                              plot_dir = paste0("./results/Personalized_analysis/20240119/boxplot/", .x, "/"),
                              plot_suffix = paste0("-", item[1], "-VS-", item[2], ".pdf"),
                              comparisons_list = list(c(item[1], item[2])),
                              x_order = c(item[1], item[2]),plot_width = 25,ncol = 6
    )
})



#### target gene ###
spatial.circ.boxplot <- function(g, PlotsPath, object, genelist){
  PlotsPath <- PlotsPath
  dir.create(PlotsPath)
  object <- exp_group_object
  num_AOI <- nrow(object@anno)
  if("SlideName" %in% names(object@anno)){
    unique_slide <- unique(object@anno$SlideName)
    num_slide <- length(unique_slide)
  }
  # 循环title
  template <- "#### 差异分组: %s {.unlisted .unnumbered}
"
  title_template <- "### 玻片: %s {.tabset .tabset-fade .unnumbered}
"
  for (per_slide in unique_slide){
    cat(sprintf(title_template, paste0(per_slide)))
    
    # 内圈循环，以差异分析为单位
    
    # 提取差异分析中最显著的基因个数
    sig_gene_list <- genelist
    # 提取本玻片上的AOI队列
    AOI_list <- rownames(subset(object@anno, SlideName == per_slide))
    
    # 开始绘图模块
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(plotrix))
    suppressPackageStartupMessages(library(ggforce))
    
    # 提取坐标
    coordinates <- t(object@anno[,c("ROICoordinateX", "ROICoordinateY")])
    # colnames(coordinates) <- colnames(coordinates) %>% str_extract("(?<=WTA\\_).*")
    # 筛选矩阵只保留显著10个基因
    exp <- object@exp[rownames(object@exp) %in% sig_gene_list, ]
    # colnames(exp) <- colnames(exp) %>% str_extract("(?<=WTA\\_).*")
    plotting_data <- as.data.frame(t(rbind(coordinates, exp)))
    # 过滤样本只保留本玻片的样本
    plotting_data <- plotting_data[rownames(plotting_data) %in% AOI_list, ]
    plotting_data$AOI <- rownames(plotting_data)
    plotting_data <- plotting_data %>% tidyr::pivot_longer(-c(ROICoordinateX,ROICoordinateY, AOI))
    plotting_data <- plotting_data %>% group_by(AOI) %>% 
      mutate(start_location = (1:n()-1) /n()*pi*2, 
             end_location = 1:n() /n()*pi*2)
    check_overplot <- as.matrix(unique(data.frame(
      ROICoordinateX = plotting_data$ROICoordinateX,
      ROICoordinateY = plotting_data$ROICoordinateY)))
    class(check_overplot) <- "numeric"
    plotting_data <- data.frame(check_overplot) %>% full_join(plotting_data)
    plotting_data$AOI <- plotting_data$AOI %>% str_extract("(?<=WTA\\_).*")
    p <- ggplot(plotting_data) + geom_arc_bar(aes(x0 = ROICoordinateX, y0 = -ROICoordinateY,
                                                  r0 = 0, r = sqrt(value * 8e7 / max(value)), start = start_location,
                                                  end   = end_location, fill  = name), alpha = 0.65, color = "transparent") + 
      # geom_text(aes(x = ROICoordinateX, y = -ROICoordinateY, label = AOI), size = 2,
      #           hjust = 0.5, vjust = 3, check_overlap = T) +
      coord_fixed() + theme_bw() +
      scale_fill_manual("Genes", values = object@config$pan_color) +
      xlab("Slide_Coordinate_X") + ylab("Slide_Coordinate_Y") +
      theme(panel.grid=element_blank())
    subchunkify(p, fig_height = 15, fig_width = 12)
    ggsave(paste0(PlotsPath, "Spatial_expression_circular_barplot_", gsub(" ", "_", g),
                  "_on_slide_", per_slide, ".pdf"), plot = p, height = 15, width = 12)
  }
}




PlotsPath <- "./results/Personalized_analysis/20240119/spatial-circ-boxplot/"
num_AOI <- nrow(exp_group_object@anno)
if("SlideName" %in% names(exp_group_object@anno)){
  unique_slide <- unique(exp_group_object@anno$SlideName)
  num_slide <- length(unique_slide)
}

# 循环title
template <- "#### 差异分组: %s {.unlisted .unnumbered}
"
title_template <- "### 玻片: %s {.tabset .tabset-fade .unnumbered}
"



# 外圈循环， 以玻片为单位
for (per_slide in unique_slide){
  cat(sprintf(title_template, paste0(per_slide)))
  
  # 内圈循环，以差异分析为单位
  for (g in names(object@deg)){
    cat(sprintf(template, paste0(object@deg[[g]]@label_string)))
    
    # 提取差异分析中最显著的基因个数
    sig_gene_list <- rownames(subset(object@deg[[g]]@result_df, 
                                     change != 'NOT') %>% slice_min(PValue, n = 10))
    # 提取本玻片上的AOI队列
    AOI_list <- rownames(subset(object@anno, SlideName == per_slide))
    
    # 筛选差异分析分组包含的AOI
    AOI_list<- AOI_list[AOI_list %in% object@deg[[g]]@filtered_sample]
    
    # 开始绘图模块
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(plotrix))
    suppressPackageStartupMessages(library(ggforce))
    
    # 提取坐标
    coordinates <- t(object@anno[,c("ROICoordinateX", "ROICoordinateY")])
    # colnames(coordinates) <- colnames(coordinates) %>% str_extract("(?<=WTA\\_).*")
    # 筛选矩阵只保留显著10个基因
    exp <- object@exp[rownames(object@exp) %in% sig_gene_list, ]
    # colnames(exp) <- colnames(exp) %>% str_extract("(?<=WTA\\_).*")
    plotting_data <- as.data.frame(t(rbind(coordinates, exp)))
    # 过滤样本只保留本玻片的样本
    plotting_data <- plotting_data[rownames(plotting_data) %in% AOI_list, ]
    plotting_data$AOI <- rownames(plotting_data)
    plotting_data <- plotting_data %>% tidyr::pivot_longer(-c(ROICoordinateX,ROICoordinateY, AOI))
    plotting_data <- plotting_data %>% group_by(AOI) %>% 
      mutate(start_location = (1:n()-1) /n()*pi*2, 
             end_location = 1:n() /n()*pi*2)
    check_overplot <- as.matrix(unique(data.frame(
      ROICoordinateX = plotting_data$ROICoordinateX,
      ROICoordinateY = plotting_data$ROICoordinateY)))
    class(check_overplot) <- "numeric"
    plotting_data <- data.frame(check_overplot) %>% full_join(plotting_data)
    plotting_data$AOI <- plotting_data$AOI %>% str_extract("(?<=WTA\\_).*")
    p <- ggplot(plotting_data) + geom_arc_bar(aes(x0 = ROICoordinateX, y0 = -ROICoordinateY,
                                                  r0 = 0, r = sqrt(value * 8e7 / max(value)), start = start_location,
                                                  end   = end_location, fill  = name), alpha = 0.65, color = "transparent") + 
      # geom_text(aes(x = ROICoordinateX, y = -ROICoordinateY, label = AOI), size = 2,
      #           hjust = 0.5, vjust = 3, check_overlap = T) +
      coord_fixed() + theme_bw() +
      scale_fill_manual("Genes", values = object@config$pan_color) +
      xlab("Slide_Coordinate_X") + ylab("Slide_Coordinate_Y") +
      theme(panel.grid=element_blank())
    subchunkify(p, fig_height = 15, fig_width = 12)
    ggsave(paste0(PlotsPath, "Spatial_expression_circular_barplot_", gsub(" ", "_", g),
                  "_on_slide_", per_slide, ".pdf"), plot = p, height = 15, width = 12)
  }
}

lapply(names(listgene), \(.x){
  PlotsPath <- paste0("./results/Personalized_analysis/20240119/spatial-circ-boxplot/", .x, "/")
  spatial.circ.boxplot(g = .x, PlotsPath = PlotsPath, object = exp_group_object, genelist = listgene[[.x]])
})

