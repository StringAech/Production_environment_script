suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(rstatix))
`%!in%` <- Negate(`%in%`)

adjusted_t_test <- function(df_melt, x_aes, comparisons_list = NULL, adjust_method = "bonferroni", output_path){
  # Pairwise t-test between groups
  # browser()
  x_aes <- x_aes
  adjust_method <- adjust_method
  stat.test <- df_melt %>%
    group_by(Dimensions) %>%
    # Pairwise t-test between groups
    # If you want to assume the equality of variances (Student t-test), specify the option var.equal = TRUE:
    # https://www.datanovia.com/en/lessons/t-test-in-r/
    t_test(
      as.formula(paste("values ~", x_aes)), 
      comparisons = comparisons_list,
      detailed = TRUE) %>%
    adjust_pvalue(method = adjust_method) %>% 
    add_significance("p")
  stat.df <- data.frame(stat.test)
  # Write the dataframe to a CSV file
  write.csv(stat.df, output_path, row.names = FALSE)
  return(as.data.frame(stat.df))
}

adjusted_wilcox_test <- function(df_melt, x_aes, comparisons_list = NULL, adjust_method = "bonferroni", output_path){
  x_aes <- x_aes
  adjust_method <- adjust_method
   # Pairwise t-test between groups
  stat.test <- df_melt %>%
    group_by(Dimensions) %>%
    # Pairwise t-test between groups
    # If you want to assume the equality of variances (Student t-test), specify the option var.equal = TRUE:
    # https://www.datanovia.com/en/lessons/t-test-in-r/
    rstatix::wilcox_test(
      as.formula(paste("values ~", x_aes)),
      comparisons = comparisons_list,
      # var.equal = TRUE, 
      detailed = TRUE) %>%
    adjust_pvalue(method = adjust_method) %>% 
    add_significance("p")
  stat.df <- data.frame(stat.test)
  # Write the dataframe to a CSV file
  write.csv(stat.df, output_path, row.names = FALSE)
  return(as.data.frame(stat.df))
}

rstatix_t_test <- function(df, x_aes, values){
  # browser()
  stat.test <- df %>%
    t_test(
      as.formula(paste(values, '~', x_aes, collapse = " ")), 
      detailed = TRUE) %>%
    add_xy_position(x = x_aes)
  return(stat.test)
}
rstatix_wilcox_test <- function(df, x_aes, values){
  stat.test <- df %>%
    rstatix::wilcox_test(
      as.formula(paste(values, '~', x_aes, collapse = " ")),
      detailed = TRUE) %>%
    add_xy_position(x = x_aes)
  return(stat.test)
}

single_boxplot_rstatix <- function(df, id_colnames, x_aes, ylab = '', xlab = '', x_order = NULL, 
                                   plot_width, ncol, plot_rownum ,label = 'p', plot_dir, plot_suffix, stat_method = stat_method){
  dir.create(plot_dir, showWarnings = F, recursive = T)
  if(is.null(x_order)) x_order <- unique(df[[x_aes]])
# browser()
  exp_colnames <- setdiff(colnames(df), id_colnames)
  for (column in exp_colnames) {
    print(column[1])
    exp_df <- dplyr::select(df, all_of(c(id_colnames, column)))
    exp_df[[x_aes]] <- factor(exp_df[[x_aes]], levels = x_order)
    switch (stat_method,
            `t.test` = {t.stat = rstatix_t_test(exp_df, x_aes, column)}, 
            `wilcox.test` = { t.stat = rstatix_wilcox_test(exp_df)})
    p <- ggboxplot(data = exp_df,
                   x = x_aes,
                   y = column,
                   color = x_aes,
                   ylab = ylab,
                   xlab = xlab,
                   order = x_order) +
      stat_pvalue_manual(t.stat, label = label, tip.length = 0)
    column <- gsub(':', '', column)
    column <- gsub('%', '', column)
    ggsave(p, filename = paste0(plot_dir, column, plot_suffix), dpi = 600, limitsize = FALSE, 
           width = plot_width,
           height = plot_width / ncol * 1.5 * plot_rownum)
  }
}

adjusted_boxplot_nojitter <- function(input_df, id_colnames, x_aes, plot_dir, comparisons_list = NULL, ylab = '', x_order = NULL,
                                      xlab = '', label = 'p', plot_suffix, plot_width = 20, ncol = 6, adjust_method = "bonferroni", stat_method = "t.test"){
  if (adjust_method %!in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")){
    stop("the adjust_method not in holm, hochberg, hommel, bonferroni, BH, BY, fdr, none \n
    please check it")
  }
  if (stat_method %!in% c("t.test","wilcox.test")){
    stop("the stat_method not t.test or wilcox.test, please check it")
  }
  # browser()
  dir.create(plot_dir, showWarnings = F, recursive = T)
  df <- input_df
  dt_melt <- pivot_longer(data = df, cols = -id_colnames, values_to = 'values', names_to = 'Dimensions')
  ## 计算差异
  switch (stat_method,
          't.test' = {stat_df <- adjusted_t_test(df_melt = dt_melt,
                                                 x_aes, 
                                                 comparisons_list = comparisons_list, 
                                                 output_path = paste0(plot_dir, 'ttest.adjusted_p.csv'),
                                                 adjust_method = adjust_method)},
          'wilcox.test' = {stat_df <- adjusted_wilcox_test(df_melt = dt_melt,
                                                           x_aes,
                                                           comparisons_list = comparisons_list, 
                                                           output_path = paste0(plot_dir, 'wilcox.adjusted_p.csv'),
                                                           adjust_method = adjust_method)}
  )
  ## 筛选差异
  significant_stat_df <- dplyr::filter(stat_df, p < 0.05)
  significant_dimensions <- significant_stat_df$Dimensions %>% unique
  insignificant_stat_df <- dplyr::filter(stat_df, p >= 0.05)
  insignificant_dimensions <-  insignificant_stat_df$Dimensions %>% unique
  insignificant_dimensions <- insignificant_dimensions[insignificant_dimensions %!in% significant_dimensions]
  ## 显著画图
  ### 设置长宽
  data_colnum <- ncol(df) - 1  # Exclude the ID column
  plot_rownum <- data_colnum / ncol
  ### 
  single_boxplot_rstatix(df = dplyr::select(df, all_of(c(id_colnames, significant_dimensions))), 
                         id_colnames = id_colnames, 
                         x_aes = x_aes, ylab = ylab, xlab = xlab, x_order = x_order, stat_method = stat_method,
                         plot_width = plot_width, ncol = ncol, label = label, 
                         plot_rownum = plot_rownum, plot_dir = paste0(plot_dir,"significant/"), plot_suffix = plot_suffix)
  single_boxplot_rstatix(df = dplyr::select(df, all_of(c(id_colnames, insignificant_dimensions))), 
                         id_colnames = id_colnames, 
                         x_aes = x_aes, ylab = ylab, xlab = xlab, x_order = x_order, stat_method = stat_method,
                         plot_width = plot_width, ncol = ncol, label = label, 
                         plot_rownum = plot_rownum, plot_dir = paste0(plot_dir,"insignificant/"), plot_suffix = plot_suffix)

}

