

exp <- exp_group_object@exp %>% t
anno <- exp_group_object@anno %>% 
  dplyr::select(CancerType)

time.point <- c("LPA", "APA", "PPA", "CAA", "SPA", "MPA")

time_arrange <- function(data, timepoint, anno.col){
  gene <- colnames(data)[-1]
  data <- data %>% 
    # group_by({{anno.col}}) %>% 
    summarise(across(everything(), mean), .by = {{anno.col}})
  rise_gene <- pblapply(gene, \(.x){
    data1 <- data %>%
      dplyr::select({{anno.col}}, .x) %>%
      arrange(.x)
    if (all(data1[, anno.col] == timepoint)) return(.x)
  })
  return(Filter(Negate(is.null), rise_gene))
}

undebug(time_arrange)

test <- merge(anno, exp, by = 0) %>% 
  dplyr::select(-1) 
time.point_rev <- rev(time.point)
library(pbapply)
start_time <- Sys.time()
time_arrange(test, timepoint = time.point_rev, anno.col = "CancerType")
end_time <- Sys.time()
print(end_time - start_time)
