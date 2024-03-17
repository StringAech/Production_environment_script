#### check the gene expression in the diff group and plot the boxpolt ####
Gene.expr.check <- function(exp_group_object, g_complex = 1, gene = "CD68"){
  deg_object_list <- setup_deg_object_from_string(exp_group_object@deg_groups)
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
  groupinfo <- forcats::fct_relevel(groupinfo,c(.control,.case))
  count <- deg_object_list[[g_complex]]@filtered_matrix
  return(boxplot(as.numeric(count[gene, ])~groupinfo))
}
Gene.expr.check(exp_group_object, 9, "TGFB1")


#### checck the gene in all group information ####
lapply(1:9, \(.x){
  exp_group_object@deg[[.x]]@result_df %>% dplyr::filter(rownames(.) == "TGFB1") %>% 
    mutate(Group = names(exp_group_object@deg)[[.x]])
}) %>% do.call(rbind, .) %>%
  rownames_to_column("Gene") %>% 
  mutate(Gene = "TGFB1")


