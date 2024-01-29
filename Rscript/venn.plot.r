library(dplyr)
library(tidyverse)
library(VennDiagram)
#### plasma ####
load("~/work/Olink/YKKY0100/Olink_plasma/results_plasma/rendering_cache_Olink/DEG-calculation_55d5c4ba314cbbb7e3ee0dca4a9def16.RData")
deg <- exp_group_object@deg %>% 
  lapply(\(.x){
    .data <- .x@result_df
    .data <- rownames_to_column(.data = .data, "Gene") %>% 
      dplyr::filter(Change != "NOT")
  }) %>% lapply(\(.x) return(if(length(.x[, 1]>0)) .x[, 1]))
deg <- deg[! sapply(deg, is.null)]# del the group which number of DEG is zero
# venn.diagram(deg[1:5], filename = NULL, category.names = names(deg)[1:5]) %>% grid.draw()
combn_list <- combn(names(deg), 2, simplify = F)
lapply(combn_list, function(.x){
  p <- venn.diagram(list(deg[[.x[1]]], deg[[.x[2]]]), 
               category.names = c(.x[1], .x[2]), imagetype = "tiff",
               alpha = .8, cat.pos = -10, cat.default.pos = "text",
               filename = NULL,
               # filename = paste0("./test/tiff/", .x[1], "_", .x[2], "venn.plot.tiff"),
               fill = col_vector[2:3], compression = "lzw",
               disable.logging = F)
  pdf(paste0("./Olink_plasma/venn.plot/", .x[1], "_", .x[2], "venn.plot.pdf"))
  grid.draw(p)
  dev.off()
})
#### CSF ####
load("~/work/Olink/YKKY0100/Olink_CSF/results_CSF/rendering_cache_Olink/DEG-calculation_721c515ef58d91962b9a95096dd3b72b.RData")
deg <- exp_group_object@deg %>% 
  lapply(\(.x){
    .data <- .x@result_df
    .data <- rownames_to_column(.data = .data, "Gene") %>% 
      dplyr::filter(Change != "NOT")
  }) %>% lapply(\(.x) return(if(length(.x[, 1]>0)) .x[, 1]))
deg <- deg[! sapply(deg, is.null)]# del the group which number of DEG is zero
# venn.diagram(deg[1:5], filename = NULL, category.names = names(deg)[1:5]) %>% grid.draw()
combn_list <- combn(names(deg), 2, simplify = F)
lapply(combn_list, function(.x){
  p <- venn.diagram(list(deg[[.x[1]]], deg[[.x[2]]]), 
                    category.names = c(.x[1], .x[2]), imagetype = "tiff",
                    alpha = .8, cat.pos = -10, cat.default.pos = "text",
                    filename = NULL,
                    # filename = paste0("./test/tiff/", .x[1], "_", .x[2], "venn.plot.tiff"),
                    fill = col_vector[2:3], compression = "lzw",
                    disable.logging = F)
  pdf(paste0("./Olink_CSF/venn.plot/", .x[1], "_", .x[2], "venn.plot.pdf"))
  grid.draw(p)
  dev.off()
})