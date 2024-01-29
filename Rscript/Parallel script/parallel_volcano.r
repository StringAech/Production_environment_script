# plot_vol_function <- function(g){
#     p <- suppressMessages(volcano_edgeR_plot(exp_group_object@deg[[g]]@result_df, genelist = NULL, 
#                                              FDRcutoff = exp_group_object@config$deg$fdr_cutoff,
#                                              lfc = exp_group_object@config$deg$lfc_cutoff,
#                                              num_label_each = exp_group_object@config$deg$volcano_num_label_each))
#   ggsave(paste0(PlotsPath, "DEG_volcano_", gsub("的差异","",exp_group_object@deg[[g]]@label_string),".pdf"), plot = p, width = 9, height = 9)
# }
# 
# # PlotsPath <- "./results_change_logfc/new_group/volcano/"
# library(parallel);#加载并行计算包
# cl <- makeCluster(10);# 初始化cpu集群
# clusterExport(cl,'exp_group_object');
# # clusterExport(cl,'CachePath');
# # clusterExport(cl,'space');
# # clusterExport(cl,'PlotsPath');
# # clusterExport(cl,'glm_diff');
# clusterEvalQ(cl,library(dplyr));#添加并行计算中用到的包
# clusterEvalQ(cl,library(EPARS));#添加并行计算中用到的包
# get1=parLapply(cl, names(exp_group_object@deg), plot_vol_function)
# parallel::stopCluster(cl)
# 
# parallel_DEG.r


parallel_volcano <- function(exp_group_object, PlotsPath, core.num = 5){
  suppressPackageStartupMessages(library(EPARS))
  suppressPackageStartupMessages(library(ggplot2))
  plot_vol_function <- function(g){
    p <- suppressMessages(volcano_edgeR_plot(exp_group_object@deg[[g]]@result_df, 
                                             genelist = NULL, 
                                             FDRcutoff = exp_group_object@config$deg$fdr_cutoff,
                                             lfc = exp_group_object@config$deg$lfc_cutoff,
                                             num_label_each = exp_group_object@config$deg$volcano_num_label_each))
    ggsave(paste0(PlotsPath, 
                  "DEG_volcano_", gsub("的差异","",exp_group_object@deg[[g]]@label_string),".pdf"),
           plot = p, width = 9, height = 9)
  }
  library(parallel);#加载并行计算包
  cl <- makeCluster(core.num);# 初始化cpu集群
  clusterExport(cl,'exp_group_object');
  # clusterExport(cl,'PlotsPath');
  # clusterExport(cl,'space');
  # clusterExport(cl,'PlotsPath');
  # clusterExport(cl,'glm_diff');
  clusterEvalQ(cl,library(dplyr));#添加并行计算中用到的包
  clusterEvalQ(cl,library(EPARS));#添加并行计算中用到的包
  get1=parLapply(cl, names(exp_group_object@deg), plot_vol_function)
  parallel::stopCluster(cl)
}


#### eg :
# parallel_volcano(exp_group_object, PlotsPath = "./results/新增加的part4/test/volcano_lplot/")
