parallel_heatmap <- function(exp_group_object, CachePath, PlotsPath, core.num = 10, pdf_height = NULL, pdf_width = NULL, ...){
  # height <- height
  # width <- width
  heatmap_function <- function(g, width_ = pdf_width, height_ = pdf_height){
    exp_group_object_single = exp_group_object
    exp_group_object_single@exp <- exp_group_object_single@exp[ , which(colnames(exp_group_object_single@exp) %in% (exp_group_object@deg[[g]]@filtered_sample))]
    exp_group_object_single@anno <- exp_group_object_single@anno[ which(rownames(exp_group_object_single@anno) %in% (exp_group_object@deg[[g]]@filtered_sample)), ]
    exp_group_object_single@heatmap_groups =  append(exp_group_object_single@heatmap_groups, 
                                                     gsub("\\(.*","", sapply(strsplit(g, ","), "[", 1)))
    deg_str <- gsub("\\(.*","", sapply(strsplit(g, ","), "[", 1))
    exp_group_object_single@anno[[deg_str]] <- exp_group_object_single@anno[[deg_str]] %>% factor()
    
    exp_group_object_single <- epars_plot_setting(exp_group_object_single)
    exp_group_object_single <- epars_heatmap_anno_data(exp_group_object_single)
    top_anno <- exp_group_object_single@plot_setting$heatmap_anno
    filtered_list <- rownames(exp_group_object@deg[[g]]@result_df %>% 
                                dplyr::filter(exp_group_object@deg[[g]]@result_df$change != "NOT"))
    mtx <- subset(exp_group_object_single@exp, rownames(exp_group_object_single@exp)
                  %in% filtered_list)
    group_df_order = order(exp_group_object_single@plot_setting$group_df[, deg_str])
    ordered_sample <- rownames(exp_group_object_single@plot_setting$group_df)[group_df_order]
    htmp_color <- exp_group_object@config$heatmap_color
    mtx <- htmp_scale(mtx)
    # system(paste("echo", g, ">>run.txt"))
    p <- Heatmap(matrix = mtx, 
                 show_column_names = F, 
                 show_row_names = T, 
                 border = T, 
                 column_order = ordered_sample, 
                 column_split = exp_group_object_single@plot_setting$group_df[, deg_str],
                 heatmap_legend_param = list(title = "", 
                                             legend_height = unit(4, "cm")), column_names_rot = -90, 
                 row_names_gp = gpar(fontsize = 8), 
                 column_names_gp = gpar(fontsize = 8), 
                 clustering_method_rows = "ward.D2", 
                 name = "Scale(exp)", 
                 cluster_column_slices = T, 
                 cluster_rows = T, 
                 # width = ncol(mtx)*unit(4.5, "mm"),
                 # height = nrow(mtx)*unit(4, "mm"),
                 cluster_columns = F, 
                 clustering_method_columns = "ward.D2",
                 top_annotation = top_anno, 
                 col = htmp_color)
    suppressMessages(save_pdf(p, paste0(PlotsPath, "DEG_heatmap_",gsub("的差异","",exp_group_object@deg[[g]]@label_string), ".pdf"), width = width_, height = height_))
  }
  CachePath <- CachePath
  PlotsPath <- PlotsPath
  library(EPARS)
  library(parallel)  #加载并行计算包
  cl <- makeCluster(core.num);  # 初始化cpu集群
  clusterExport(cl, 'exp_group_object');
  # clusterExport(cl, 'CachePath');
  # clusterExport(cl, "space")
  # clusterExport(cl, "PlotsPath");
  # clusterExport(cl, "glm_diff")
  # clusterExport(cl, "rescale.AsIs")
  clusterEvalQ(cl, library(dplyr));  #添加并行计算中用到的包
  clusterEvalQ(cl, library(EPARS)); #添加并行计算中用到的包
  clusterEvalQ(cl, library(RColorBrewer))  #添加并行计算中用到的包
  clusterEvalQ(cl, library(R.utils))  #添加并行计算中用到的包
  clusterEvalQ(cl, library(ComplexHeatmap))  #添加并行计算中用到的包
  clusterEvalQ(cl, library(tidyHeatmap))  #添加并行计算中用到的包
  enrichment.res <- parLapply(cl, names(exp_group_object@deg), heatmap_function)
  parallel::stopCluster(cl)
}

# parallel_heatmap(exp_group_object, CachePath = "./Heatmap/", PlotsPath = "./Heatmap/", core.num = 6)

