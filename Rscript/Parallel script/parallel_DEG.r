# for example :
# exp_group_object <- parallel_DEG(exp_group_object, CachePath = "./test/", core.num = 25)

obtain_edgeR_DE <- function (matrix, groupinfo, lfc_cutoff = 1, fdr_cutoff = 0.05, 
          outdir = "data", prefix = "DEGs", cutoff_sig = "FDR") 
{
  library(edgeR)
  dge <- DGEList(counts = matrix, group = groupinfo)
  keep <- filterByExpr(matrix, group = groupinfo, min.count = 3)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge$samples$lib.size <- colSums(dge$counts)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~0 + groupinfo)
  rownames(design) <- colnames(dge)
  colnames(design) <- levels(groupinfo)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit <- glmFit(dge, design)
  fit2 <- glmLRT(fit, contrast = c(-1, 1))
  DEG = topTags(fit2, n = nrow(exp))
  DEG = as.data.frame(DEG)
  logFC_cutoff <- lfc_cutoff
  switch(cutoff_sig, 
         "P_value" = {
           k1 = (DEG$PValue < fdr_cutoff) & (DEG$logFC < -logFC_cutoff)
           k2 = (DEG$PValue < fdr_cutoff) & (DEG$logFC > logFC_cutoff)
         },
         "FDR" = {
           k1 = (DEG$FDR < fdr_cutoff) & (DEG$logFC < -logFC_cutoff)
           k2 = (DEG$FDR < fdr_cutoff) & (DEG$logFC > logFC_cutoff)
         })
  DEG$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
  DEG <- DEG %>% arrange(factor(change, levels = c("UP", "DOWN", 
                                                   "NOT")), PValue)
  lcpm <- cpm(matrix, log = TRUE)
  write.table(lcpm, paste0(outdir, "/", prefix, "_TMM_normalized_edgeR.csv"), 
              sep = ",", row.names = T, col.names = NA)
  write.table(DEG, paste0(outdir, "/", prefix, "_edgeR.csv"), 
              sep = ",", row.names = T, col.names = NA)
  return(DEG)
}
parallel_DEG <- function(exp_group_object, CachePath, core.num = 20){
  obtain_edgeR_DE <- function (matrix, groupinfo, lfc_cutoff = 1, fdr_cutoff = 0.05, 
                               outdir = "data", prefix = "DEGs", cutoff_sig = "FDR") 
  {
    library(edgeR)
    dge <- DGEList(counts = matrix, group = groupinfo)
    keep <- filterByExpr(matrix, group = groupinfo, min.count = 3)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge$samples$lib.size <- colSums(dge$counts)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~0 + groupinfo)
    rownames(design) <- colnames(dge)
    colnames(design) <- levels(groupinfo)
    dge <- estimateGLMCommonDisp(dge, design)
    dge <- estimateGLMTrendedDisp(dge, design)
    dge <- estimateGLMTagwiseDisp(dge, design)
    fit <- glmFit(dge, design)
    fit2 <- glmLRT(fit, contrast = c(-1, 1))
    DEG = topTags(fit2, n = nrow(exp))
    DEG = as.data.frame(DEG)
    logFC_cutoff <- lfc_cutoff
    switch(cutoff_sig, 
           "P_value" = {
             k1 = (DEG$PValue < fdr_cutoff) & (DEG$logFC < -logFC_cutoff)
             k2 = (DEG$PValue < fdr_cutoff) & (DEG$logFC > logFC_cutoff)
           },
           "FDR" = {
             k1 = (DEG$FDR < fdr_cutoff) & (DEG$logFC < -logFC_cutoff)
             k2 = (DEG$FDR < fdr_cutoff) & (DEG$logFC > logFC_cutoff)
           })
    DEG$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
    DEG <- DEG %>% arrange(factor(change, levels = c("UP", "DOWN", 
                                                     "NOT")), PValue)
    lcpm <- cpm(matrix, log = TRUE)
    write.table(lcpm, paste0(outdir, "/", prefix, "_TMM_normalized_edgeR.csv"), 
                sep = ",", row.names = T, col.names = NA)
    write.table(DEG, paste0(outdir, "/", prefix, "_edgeR.csv"), 
                sep = ",", row.names = T, col.names = NA)
    return(DEG)
  }
  suppressMessages(library(forcats))
  suppressPackageStartupMessages(library(EPARS))
  # Based on the input string information for DEG groups, setup basic structure
  #  for deg_object_list before any DE calculation.
  deg_object_list <- setup_deg_object_from_string(exp_group_object@deg_groups)
  #### defferential gene analysis ####
  library(parallel);#加载并行计算包
  cl <- makeCluster(core.num);# 初始化cpu集群
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
    system(paste("echo",g_complex,">>run.txt"))
    
    #输入差异分组时,要求control在前,case在后
    # .group <- deg_object_list[[g_complex]]@groupinfo[[1]]
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
                                 fdr_cutoff = exp_group_object@config$deg$fdr_cutoff,
                                 outdir = CachePath,
                                 prefix = paste0("DE_list_" ,deg_object_list[[g_complex]]@label_string),
                                 cutoff_sig = "FDR")
      
      tmp_deg=data.frame(tmp_deg,check.names = F)
      write.table(tmp_deg, paste0(CachePath, "/", paste0("DE_list_" ,deg_object_list[[g_complex]]@label_string), "_edgeR.csv"),
                  sep = ",", row.names = T, col.names = NA)
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
    return(data.frame())
  }
  clusterExport(cl,'exp_group_object');
  # clusterExport(cl,'deg_object_list');
  # clusterExport(cl,'CachePath');
  clusterEvalQ(cl,library(dplyr));#添加并行计算中用到的包
  clusterEvalQ(cl,library(EPARS));#添加并行计算中用到的包
  get1=parLapply(cl,1:length(deg_object_list), do.data)
  for (g_complex in 1:length(deg_object_list)){
    deg_object_list[[g_complex]] <- get1[[g_complex]]
  }
  parallel::stopCluster(cl)
  
  saveRDS(deg_object_list, paste0(CachePath, "/deg_object_list.RDS"))
  for(i in deg_object_list%>%names){
    deg_object_list[[i]]
  }
  exp_group_object@deg<- deg_object_list
  # 给差异分析增加人类可读的描自然语言描述
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
  saveRDS(exp_group_object, paste0(CachePath, "/exp_group_object.RDS")) 
  return(exp_group_object)
}


# exp_group_object <- parallel_DEG(exp_group_object, CachePath = "./test/", core.num = 25)
