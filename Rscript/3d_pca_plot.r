plot3d_point <- function(data, anno, color = "balck", .fill = "red", pch.log = T, dir = "./", main = "abd", ...){
  ## calculate the type number
  num <- length(data[, anno] %>% unique())
  # if (! pch.log == T & num <= 5) {
    shapes <- c(21:25)[1:num] %>% .[as.numeric(as.factor(data[, anno]))]
    # stop("asdfasdfasdfasdfasdfasfaf")
  # }
  colors <- .fill[as.numeric(as.factor(data[, anno]))]
  pdf(dir)
  s3d <- scatterplot3d(dplyr::select(data, -anno), pch = shapes, box = T, color = "black", 
                       type = "p", main = main, 
                       angle = 45, asp = NA, lty.hide = 2,
                       scale.y = 1, bg = colors)
  legend("right", bg = "transparent", 
    legend = unique(data[, anno]), horiz = F, pch = unique(shapes),
    col = "black", pt.bg = unique(colors))
  dev.off()
  
  
}

plot3d_point(data[,-1], anno = "Acute", .fill = colors, dir = "./ad.pdf")

#### plasma ####
file.path <- list.files(path = "./results_plasma/2.Overall_Analysis/", pattern = "PCA.*csv")
file.dir <- paste0("./results_plasma/2.Overall_Analysis/", file.path)
dir.create("./results_plasma/2.Overall_Analysis/3d_plot")
data <- list()
for (i in file.dir){
  data[[basename(i)]] <- read.csv(i) %>% 
    dplyr::select(str_extract(basename(i), "(?<=PCA\\_).*(?=\\.csv)"), 
                  PC1, PC2, PC3)
}
# lapply(data, plot3d_point, anno = colnames(x)[1], .fill = colors, dir = "./results_plasma/2.Overall_Analysis/3d_plot/")
lapply(names(data), function(.x){
  .data <- data[[.x]]
  plot3d_point(data = .data, anno = colnames(.data)[1], .fill = colors, 
               dir = paste0("./results_plasma/2.Overall_Analysis/3d_plot/",
                            gsub(x = .x, pattern = "\\.csv", replacement = ""),
                            "_plot.pdf"
                            ),
               main = gsub(x = .x, pattern = "\\.csv", replacement = "")
               )
})
#### CSF ####
setwd("~/work/Olink/YKKY0100/Olink_CSF")
file.path <- list.files(path = "./results_CSF/2.Overall_Analysis/", pattern = "PCA.*csv")
file.dir <- paste0("./results_CSF/2.Overall_Analysis/", file.path)
dir.create("./results_CSF/2.Overall_Analysis/3d_plot")
data <- list()
for (i in file.dir){
  data[[basename(i)]] <- read.csv(i) %>% 
    dplyr::select(str_extract(basename(i), "(?<=PCA\\_).*(?=\\.csv)"), 
                  PC1, PC2, PC3)
}
# lapply(data, plot3d_point, anno = colnames(x)[1], .fill = colors, dir = "./results_plasma/2.Overall_Analysis/3d_plot/")
lapply(names(data), function(.x){
  .data <- data[[.x]]
  plot3d_point(data = .data, anno = colnames(.data)[1], .fill = colors, 
               dir = paste0("./results_CSF/2.Overall_Analysis/3d_plot/",
                            gsub(x = .x, pattern = "\\.csv", replacement = ""),
                            "_plot.pdf"
               ),
               main = gsub(x = .x, pattern = "\\.csv", replacement = "")
  )
})
