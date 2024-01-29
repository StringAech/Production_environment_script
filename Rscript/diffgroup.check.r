diff_sample_check <- function(deg_groups){
  deg_object_list <- setup_deg_object_from_string(deg_groups)
  test <- list()
  for (g_complex in 1:length(deg_object_list)){
    # print(g_complex)
    # 使用整理好的过滤信息对AOI进行过滤
    count = exp_group_object@exp
    groupinfo <- factor(data.frame(exp_group_object@anno)[[deg_object_list[[g_complex]]@groupinfo[[1]] ]])
    ROI=exp_group_object@anno$SlideName
    
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
      ROI=ROI[tmp_filtered_list]
      groupinfo <- factor(groupinfo[tmp_filtered_list])
    }
    #输入差异分组时,要求control在前,case在后
    .group <- deg_object_list[[g_complex]]@groupinfo[[1]]
    .case <- deg_object_list[[g_complex]]@groupinfo[[3]]
    .control <- deg_object_list[[g_complex]]@groupinfo[[2]]
    #按照输入的分组顺序进行排序
    groupinfo <- fct_relevel(groupinfo,c(.control,.case))
    test[[names(deg_object_list[g_complex])]] <- groupinfo
  }
  return(test)
}
suppressMessages(library(forcats))
# deg_groups <- c("Subpopulation(tumor:stroma),Treatment$naive; Subpopulation(tumor:Macrophage),Treatment$naive; Subpopulation(tumor:T_cell),Treatment$naive; Subpopulation(tumor:Tertiary_lymphatic_structure),Treatment$naive; Subpopulation(tumor:NTNM),Treatment$naive; Spatial_location(Peritumor:Intratumor),Treatment$naive&Subpopulation$stroma; Spatial_location(Peritumor:Intratumor),Treatment$naive&Subpopulation$Macrophage; Spatial_location(Peritumor:Intratumor),Treatment$naive&Subpopulation$T_cell; Spatial_location(Peritumor:Intratumor),Treatment$naive&Subpopulation$NTNM; Histological_subtype(Acinar:papillary),Treatment$naive; Histological_subtype(Acinar:micropapillary),Treatment$naive; Histological_subtype(Acinar:solid),Treatment$naive; Histological_subtype(Acinar:papillary),Subpopulation$stroma&Treatment$naive; Histological_subtype(Acinar:micropapillary),Subpopulation$stroma&Treatment$naive; Histological_subtype(Acinar:solid),Subpopulation$stroma&Treatment$naive; Histological_subtype(Acinar:papillary),Subpopulation$T_cell&Treatment$naive; Histological_subtype(Acinar:micropapillary),Subpopulation$T_cell&Treatment$naive; Histological_subtype(Acinar:solid),Subpopulation$T_cell&Treatment$naive; Histological_subtype(Acinar:papillary),Subpopulation$Macrophage&Treatment$naive; Histological_subtype(Acinar:micropapillary),Subpopulation$Macrophage&Treatment$naive; Histological_subtype(Acinar:solid),Subpopulation$Macrophage&Treatment$naive; Histological_subtype(Acinar:papillary),Subpopulation$NTNM&Treatment$naive; Histological_subtype(Acinar:micropapillary),Subpopulation$NTNM&Treatment$naive; Histological_subtype(Acinar:solid),Subpopulation$NTNM&Treatment$naive; Spatial_location(Peritumor:Intratumor),Histological_subtype$Acinar&Subpopulation$stroma&Treatment$naive; Spatial_location(Peritumor:Intratumor),Histological_subtype$papillary&Subpopulation$stroma&Treatment$naive; Spatial_location(Peritumor:Intratumor),Histological_subtype$solid&Subpopulation$stroma&Treatment$naive; Spatial_location(Peritumor:Intratumor),Histological_subtype$Acinar&Subpopulation$T_cell&Treatment$naive; Spatial_location(Peritumor:Intratumor),Histological_subtype$Acinar&Subpopulation$NTNM&Treatment$naive; sample_ID(naive1:naive1_lym),Subpopulation$tumor; sample_ID(naive1:naive1_lym),Subpopulation$stroma; sample_ID(naive1:naive1_lym),Subpopulation$T_cell; sample_ID(naive1:naive1_lym),Subpopulation$NTNM; sample_ID(naive1:naive1_lym),Subpopulation$Macrophage; Histological_subtype(Acinar:micropapillary),sample_ID$naive1&Spatial_location2$Peritumor_Intratumor; Histological_subtype(Acinar:solid),sample_ID$naive1&Spatial_location2$Peritumor_Intratumor; Histological_subtype(Acinar:papillary),sample_ID$naive2&Spatial_location2$Peritumor_Intratumor")
test <- diff_sample_check(deg_groups = exp_group_object@deg_groups)
sample.num <- lapply(test, \(.x){
  res <- list(sample.num = length(.x), level.num = length(levels(.x)))
  return(res)
})
# check level number 
# sapply(sample.num, \(.x) .x$level.num == 2) %>% all() # check level number 
if (!sapply(sample.num, \(.x) .x$level.num == 2) %>% all()){
  warning("Levels number error!")
         sapply(sample.num, \(.x) .x$level.num != 2 ) %>% sample.num[.] %>% names() %>% 
           gsub(pattern = "=", replacement = "$", x = .) %>% print()
  # sapply(sample.num, \(.x) .x$sample.num < 3) %>% sample.num[.] %>% names()
}
# check sample number 
if (!sapply(sample.num, \(.x) .x$sample.num > 2) %>% all()) {
  warning("Sample number error!")
  sapply(sample.num, \(.x) .x$sample.num < 3) %>% sample.num[.] %>% 
    names() %>% gsub(pattern = "=", replacement = "$", x = .) %>% print()
}



  