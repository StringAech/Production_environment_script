# source("~/work/DSP/CTA/YKKY0136/DSP_CTA/resource/.Rprofile")
source("~/script/Rscript/all.function.R")
library(tidyverse)
res <- readRDS("/home/rstudio/user/guohao/work/DSP/CTA/YKKY0136/20240127_之前的全部作废/DSP_CTA_1/results/Spatialdecon/spatialdecon.RDs")
deg_list <- readRDS("/home/rstudio/user/guohao/work/DSP/CTA/YKKY0136/20240127_之前的全部作废/DSP_CTA_1/deg_object_list.RDS")
exp_group_object <- readRDS("/home/rstudio/user/guohao/work/DSP/CTA/YKKY0136/20240127_之前的全部作废/DSP_CTA_1/exp_group_object.RDS")
cell_exp <- res$beta %>% t %>% as.data.frame()
cell_type <- colnames(cell_exp)
plot_funciton <- function(.x, cell_exp, dir){
  source("~/script/Rscript/all.function.R")
  message(paste0("Tidy **", .x@label_string, "** ing..."))
  # system(paste("echo", .x@label_string ,">>run.txt"))
  cell_type <- colnames(cell_exp)
  cell_exp <- merge(exp_group_object@anno, cell_exp, by = 0) %>% 
    dplyr::filter(Row.names %in% .x@filtered_sample) %>% 
    dplyr::select(Row.names, .x@groupinfo[[1]], cell_type)
  cell_exp <- cbind(cell_exp[, 1:2],
                    cell_exp[, 3:ncol(cell_exp)] %>% .[, apply(., 2, sum)!=0]
  )
  del_cell_type <- setdiff(cell_type, colnames(cell_exp)[-c(1:2)])
  id_colnames <- colnames(cell_exp)[1:2]
  # system(paste("echo",.x@filtered_sample,">>run.txt"))
  tryCatch({adjusted_boxplot_nojitter(
    cell_exp,
    id_colnames,
    id_colnames[2],
    paste0("./", dir, "/", gsub("的差异", "", .x@label_string), "/"),
    'cell_abundance_diff.boxplot.pdf',
    # list(c("Sensitive-pre", "Sensitive-post"), c("Resistant-pre", "Resistant-post"), 
    #      c("Resistant-pre", "Sensitive-pre"), c("Resistant-post", "Sensitive-post")),
    c(.x@groupinfo[[2]], .x@groupinfo[[3]]), 
    x_order = c(.x@groupinfo[[3]], .x@groupinfo[[2]]),
    # x_order = c('Resistant-post', 'Resistant-pre', 'Sensitive-post', 'Sensitive-pre'),
    25,  # plot_width
    6,  # plot ncol
  )
    plot_data <- pivot_longer(cell_exp, cols = -id_colnames, 
                              names_to = "Cell_type", values_to = "Abundance")
    plot_data[[.x@groupinfo[[1]]]] <- factor(plot_data[[.x@groupinfo[[1]]]], 
                                             levels = c(.x@groupinfo[[3]], .x@groupinfo[[2]]))
    p <- ggboxplot(data = plot_data, x = "Cell_type", 
                   y = "Abundance", ylab = "Abundance", 
                   xlab = "", fill = .x@groupinfo[[1]],
    ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
    ggsave(plot = p, filename = paste0("./", dir, "/", gsub("的差异", "", .x@label_string), "/", "all.celltype.boxplot.pdf"))
    
  },error = function(e){message(paste0("ERROR: **", .x@label_string, "** "))}, 
  finally = {
    write.xlsx(cell_exp,
               file = paste0("./",dir,"/",gsub("的差异", "", .x@label_string),"/","cell_type_abundance.xlsx"))
    if (length(del_cell_type) != 0) {
      write.table(del_cell_type, 
                  file = paste0("./", dir, "/", gsub("的差异", "", .x@label_string), "/", "del_cell_type.txt"), 
                  sep = "\t",quote = F,
                  row.names = F)
    }
  }
  )
}
# lapply(deg_list, plot_funciton, cell_exp = cell_exp, dir = "cell_type_plot")
library(parallel);#加载并行计算包
cl <- makeCluster(10);# 初始化cpu集群
clusterExport(cl,'exp_group_object');
clusterExport(cl,'deg_list');
clusterEvalQ(cl,library(dplyr));#添加并行计算中用到的包
clusterEvalQ(cl,library(EPARS));#添加并行计算中用到的包
get1=parLapply(cl, deg_list, plot_funciton, cell_exp = cell_exp, dir = "./results/cell_type_abundance_plot")
parallel::stopCluster(cl)



# for (i in names(deg_list)) {
#   plot_funciton(deg_list[[i]], cell_exp = cell_exp, dir = "./results/cell_type_abundance_plot")
# }
# 
# plot_funciton(deg_list[[13]], cell_exp = cell_exp, dir = "./results/cell_type_abundance_plot")


