# param 'Object': S4object exp_group_object
# param 'ident_col': cluster colnumn name, just like Seurat function FindAllMarker
# param 'CachePath': output path
# param 'core.num': parallel number
DSP_FindAllMarker <- function(object, ident_col, CachePath, core.num = 6, ...){
  ## set the deg_groups
  cluster <- object@anno[, ident_col] %>% unique()
  char <- lapply(cluster, \(.x){
    paste0("Group_", .x, "(Other_cluster:", .x, ")")
  })
  object@deg_groups <- do.call(paste, c(char, sep = "; "))
  ## set the anno data.frame
  for(i in cluster){
    col_name <- paste0("Group_", i)
    object@anno <- object@anno %>% 
      dplyr::mutate(!!sym(col_name) := ifelse(CancerType != {{i}}, "Other_cluster", {{i}}))
  }
  ## calculation the different expression genes
  ### to use parallel procession, source R script
  source("~/script/Rscript/Parallel script/library_EPARS_parallel.r")
  object <- parallel_DEG(object, CachePath = CachePath, core.num = core.num)
  return(object)
}

# exp_group_object_test <- DSP_FindAllMarker(exp_group_object, 
#                                            ident_col = "CancerType", 
#                                            CachePath = "./personal_analysis/rmd_Results/test/")