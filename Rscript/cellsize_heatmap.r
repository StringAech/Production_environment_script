library(stringr)
for (g in names(exp_group_object@deg)){
  
  exp_group_object_single = exp_group_object
  exp_group_object_single@heatmap_groups <- str_extract(g, ".*(?=\\()")
  exp_group_object_single@exp <- exp_group_object_single@exp[ , which(colnames(exp_group_object_single@exp) %in% (exp_group_object@deg[[g]]@filtered_sample))]
  exp_group_object_single@anno <- exp_group_object_single@anno[ which(rownames(exp_group_object_single@anno) %in% (exp_group_object@deg[[g]]@filtered_sample)), ]
  # exp_group_object_single@heatmap_groups =  gsub("\\(.*","", sapply(strsplit(g, ","), "[", 1))
  
  exp_group_object_single <- epars_plot_setting(exp_group_object_single)
  exp_group_object_single <- epars_heatmap_anno_data(exp_group_object_single)
  top_anno <- exp_group_object_single@plot_setting$heatmap_anno
  filtered_list <- rownames(exp_group_object@deg[[g]]@result_df %>% 
                              dplyr::filter(exp_group_object@deg[[g]]@result_df$change != "NOT"))
  mtx <- subset(exp_group_object_single@exp, rownames(exp_group_object_single@exp)
                %in% filtered_list)
  mtx <- htmp_scale(mtx)
  htmp_color <- exp_group_object_single@config$heatmap_color
  # group_df_order = order(exp_group_object_single@plot_setting$group_df[,2])
  # ordered_sample <- rownames(exp_group_object_single@plot_setting$group_df)[group_df_order]
  p <- Heatmap(matrix = mtx, 
               show_column_names = T,  
               show_row_names = T, border = T, 
               # column_order = ordered_sample, 
               # column_split = exp_group_object_single@plot_setting$group_df[,1],
               heatmap_legend_param = list(title = "", 
                                           legend_height = unit(4, "cm")), column_names_rot = -90, 
               row_names_gp = gpar(fontsize = 6), 
               column_names_gp = gpar(fontsize = 8), 
               clustering_method_rows = "ward.D2", 
               name = "Scale(exp)", 
               cluster_column_slices = T, 
               cluster_rows = T, 
               cluster_columns = T, 
               width = ncol(mtx)*unit(2.5, "mm"),
               height = nrow(mtx)*unit(2, "mm"),
               clustering_method_columns = "ward.D2",
               top_annotation = top_anno, 
               col = htmp_color)
  
save_pdf(p, 
         paste0("~/work/DSP/YKKY0153广东省人民医院梅老师/20231220_change_the_group/DSP_WTA/results/new_heatmap/DEG_heatmap_",gsub("的差异","",exp_group_object@deg[[g]]@label_string), ".pdf"), 
         width =  10, 
         height = 14)

}
