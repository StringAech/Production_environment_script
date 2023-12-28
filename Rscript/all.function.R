# load libraries
library(readxl)
library(rstudioapi)
library(tidyr)
library(broom)
library(dplyr)
#library(writexls)
library(heatmap3)
#library(xlsx)
library(openxlsx)
library(ComplexHeatmap)
library(reshape2)
library(RColorBrewer)
library(corrplot)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(plyr)
library(data.table)
library(hash)
library(psych)
library(readr)
library(ggsci)
library(ggsignif)
library(ggpubr)  # ggadjust_pvalue.R
library(rstatix)  # multiple group test
# Plot paired data
# library(PairedData)
`%!in%` <- Negate(`%in%`)
# create dir if not exist
make_dir <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

# convert data.frame to numeric
numeric_convert <- function(df){
  df[,1:ncol(df)] = lapply(
    1:ncol(df), function(x) {
      tryCatch({
        as.numeric(df[[x]])
      },warning = function(w) {
        df[[x]]}
      )} 
  )
  # Fill NaN values with zero
  df[is.na(df)] <- 0
  return(df)
}

numeric_only_df <- function(df){
  # Assuming your dataframe is called 'df'
  numeric_columns <- df %>%
    select_if(is.numeric)
  numeric_df <- as.data.frame(numeric_columns)
  return(numeric_df)
}

filter_df <- function(my_cols, df){
  # Filter the data frame to keep only the specified columns
  filtered_df <- subset(df, select = my_cols)
  return(filtered_df)
}

filter_row <- function(key_column, keyword, df){
  # Filter the data frame rows according to column value
  # filtered_df <- df[df$group1 == keyword,]
  filtered_df <- filter(df, grepl(keyword, df[[key_column]], ignore.case = FALSE))
  return(filtered_df)
}

filter_list <- function(keyword, input_list){
  # Filter the list to keep only the strings containing the interested keyword
  filtered_list <- input_list[grep(keyword, input_list)]
  return(filtered_list)
}

boxplot <- function(df, id_colnames, x_aes, plot_path, comparisons_list, plot_width, ncol){
  suppressWarnings({
    dt_melt <- melt(data.table(df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
    data_colnum <- ncol(df) - 1  # Exclude the ID column
    plot_rownum <- data_colnum / ncol
    p <- ggplot(data = dt_melt,
                # aes_string(my_aes)
                aes(
                  x=.data[[x_aes]], 
                  y=values, 
                  color=.data[[x_aes]]
                )
    ) + geom_boxplot() + 
      theme_classic2() +
      theme(legend.position = "none") +  # no color legend 
      #scale_fill_jco() + 
      facet_wrap(~ Dimensions, scales="free", ncol = ncol) +
      # facet_grid(~ Dimensions, scales="free", , space='free') +
      stat_compare_means(aes(method = "t.test", label = "p = {p.format}"), comparisons = comparisons_list) +
      geom_jitter(shape=16, position=position_jitter(0.2))
    ggsave(
      plot_path, 
      width = plot_width, 
      height = plot_width / ncol * 1.5 * plot_rownum,  # as.integer(float_value)
      dpi = 600,
      limitsize = FALSE
    )
  })
}
boxplot_nojitter <- function(df, id_colnames, x_aes, plot_path, comparisons_list, plot_width, ncol){
  suppressWarnings({
    dt_melt <- melt(data.table(df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
    
    data_colnum <- ncol(df) - 1  # Exclude the ID column
    plot_rownum <- data_colnum / ncol
    p <- ggplot(data = dt_melt,
                # aes_string(my_aes)
                aes(
                  x=.data[[x_aes]], 
                  y=values, 
                  color=.data[[x_aes]]
                )
    ) + geom_boxplot() + 
      theme_classic2() +
      theme(legend.position = "none") +  # no color legend 
      #scale_fill_jco() + 
      facet_wrap(~ Dimensions, scales="free", ncol = ncol) +
      # facet_grid(~ Dimensions, scales="free", , space='free') +
      stat_compare_means(aes(method = "t.test", label = "p = {p.format}"), comparisons = comparisons_list) + 
      stat_compare_means() # Add global Anova p-value
    ggsave(
      plot_path, 
      width = plot_width, 
      height = plot_width / ncol * 1.5 * plot_rownum,  # as.integer(float_value)
      dpi = 600,
      limitsize = FALSE
    )
  })
}
adjusted_boxplot_nojitter <- function(input_df, id_colnames, x_aes, plot_dir, 
                                      plot_suffix, comparisons_list, x_order, 
                                      plot_width, ncol, stat_method = "wilcox.test", 
                                      adjust_method = "bonferroni"
){
  if (adjust_method %!in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")){
    stop("the adjust_method not in holm, hochberg, hommel, bonferroni, BH, BY, fdr, none \n
    please check it")
  }
  if (stat_method %!in% c("t.test","wilcox.test")){
    stop("the stat_method not t.test or wilcox.test, please check it")
  }
  make_dir(plot_dir)
  df <- input_df
  names(df)[names(df) == x_aes] <- "group"
  # Replace "banana" with "grape"
  id_colnames[id_colnames == x_aes] <- "group"
  plot_path <- paste0(plot_dir, plot_suffix)
  dt_melt <- melt(data.table(df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
  # nomelt_boxplot_nojitter(df, dt_melt, id_colnames, "group", plot_path, comparisons_list, plot_width, ncol)
  
  switch (stat_method,
          `t.test` = {stat_df <- adjusted_t_test(dt_melt, id_colnames, x_aes, paste0(plot_dir, 'ttest.adjusted_p.csv'), adjust_method = adjust_method)}, 
          `wilcox.test` = {stat_df <- adjusted_wilcox_test(dt_melt, id_colnames, x_aes, paste0(plot_dir, 'wilcox.adjusted_p.csv'), adjust_method = adjust_method)}
  )
  
  # stat_df <- adjusted_wilcox_test(dt_melt, id_colnames, x_aes, paste0(plot_dir, 'adjusted_p.csv'))
  
  significant_stat_df <- subset(stat_df, p < 0.05)
  significant_stat_df$Dimensions <- as.character(significant_stat_df$Dimensions)
  significant_dimensions <- as.vector(unique(significant_stat_df$Dimensions))
  
  insignificant_stat_df <- subset(stat_df, p >= 0.05)
  insignificant_stat_df$Dimensions <- as.character(insignificant_stat_df$Dimensions)
  insignificant_dimensions <- as.vector(unique(insignificant_stat_df$Dimensions))
  insignificant_dimensions <- insignificant_dimensions[insignificant_dimensions %!in% significant_dimensions]
  single_boxplot_rstatix(
    filter_df(c(id_colnames, significant_dimensions), df),
    id_colnames, "group", paste0(plot_dir, 'significant/'), plot_suffix, x_order, stat_method = stat_method
  )
  single_boxplot_rstatix(
    filter_df(c(id_colnames, insignificant_dimensions), df),
    id_colnames, "group", paste0(plot_dir, 'insignificant/'), plot_suffix, x_order, stat_method = stat_method
  )
  # single_boxplot_nojitter(
  #     filter_df(c("group", significant_dimensions), df), 
  #     id_colnames, "group", plot_path, comparisons_list)
  # output_path <- sub('boxplot.pdf', 'paired.csv', plot_path)
  # adjusted_paired_wilcox_test(df, id_colnames, x_aes, output_path, comparisons_list)
}
single_boxplot_rstatix <- function(df, id_colnames, x_aes, plot_dir, plot_suffix, x_order,stat_method){
  make_dir(plot_dir)
  suppressWarnings({
    exp_colnames <- setdiff(colnames(df), id_colnames)
    for (column in exp_colnames) {
      print(column[1])
      exp_df <- as.data.frame(filter_df(c(id_colnames, column), df))
      names(exp_df)[names(exp_df) == column] <- "exp"
      exp_df$group <- factor(exp_df$group,levels = x_order)
      switch (stat_method,
              `t.test` = {t.stat = rstatix_t_test(exp_df)}, 
              `wilcox.test` = { t.stat = rstatix_wilcox_test(exp_df)}
      )
      # t.stat = rstatix_wilcox_test(
      #   # melt(data.table(exp_df), variable.name="group", value.name = "exp")
      #   exp_df
      # )
      p <- ggboxplot(data = exp_df,
                     x="group",
                     y="exp",
                     color="group",
                     ylab = "values",
                     xlab = "group",
                     order = x_order
      ) +
        stat_pvalue_manual(t.stat, label = "p", tip.length = 0)
      #  labs(subtitle = get_test_label(t.stat, detailed = TRUE))
      column <- gsub(':', '', column)
      column <- gsub('%', '', column)
      ggsave(
        paste0(plot_dir, column, plot_suffix), 
        # width = plot_width, 
        # height = plot_width / ncol * 1.5 * plot_rownum,  # as.integer(float_value)
        dpi = 600,
        limitsize = FALSE
      )
    }
  })
}

# ggboxplot with default t.test as method for calcualting p values
single_boxplot_nojitter <- function(df, id_colnames, x_aes, plot_path, comparisons_list){
  suppressWarnings({
    exp_colnames <- setdiff(colnames(df), id_colnames)
    for (column in exp_colnames) {
      print(column[1])
      exp_df <- as.data.frame(filter_df(c(x_aes, column), df))
      names(exp_df)[names(exp_df) == column] <- "exp"
      p <- ggboxplot(data = exp_df,
                     x="group",
                     y="exp",
                     color="group"
      ) + 
        # facet_wrap(~ Dimensions, scales="free", ncol = ncol) +
        stat_compare_means(
          aes(method = "t.test", label = "p = {p.format}"), 
          comparisons = comparisons_list)
      column <- gsub(':', paste0(''), column)
      ggsave(
        gsub('.boxplot.pdf', paste0('.', column, '.boxplot.pdf'), plot_path), 
        # width = plot_width, 
        # height = plot_width / ncol * 1.5 * plot_rownum,  # as.integer(float_value)
        dpi = 600,
        limitsize = FALSE
      )
    }
  })
}

nomelt_boxplot_nojitter <- function(df, dt_melt, id_colnames, x_aes, plot_path, comparisons_list, plot_width, ncol){
  suppressWarnings({
    data_colnum <- ncol(df) - 1  # Exclude the ID column
    plot_rownum <- data_colnum / ncol
    p <- ggboxplot(data = dt_melt,
                   x=x_aes, 
                   y='values', 
                   color=x_aes
    ) + 
      facet_wrap(~ Dimensions, scales="free", ncol = ncol) +
      # facet_grid(~ Dimensions, scales="free", , space='free') +
      stat_compare_means(
        aes(method = "t.test", label = "p = {p.format}"), 
        comparisons = comparisons_list) + 
      stat_compare_means() # Add global Anova p-value
    ggsave(
      plot_path, 
      width = plot_width, 
      height = plot_width / ncol * 1.5 * plot_rownum,  # as.integer(float_value)
      dpi = 600,
      limitsize = FALSE
    )
  })
}

paired_rstatix_t_test <- function(df) {
  stat.test <- df %>%
    # group_by(group) %>%
    t_test(
      exp ~ group, 
      # var.equal = TRUE, 
      detailed = TRUE,
      paired = TRUE
    ) %>%
    add_xy_position(x = "group")%>% 
    adjust_pvalue() %>% 
    add_significance("p.adj")
  return(stat.test)
}

paired_rstatix_wilcox_test <- function(df) {
  stat.test <- df %>%
    # group_by(group) %>%
    rstatix::wilcox_test(
      exp ~ group, 
      # var.equal = TRUE, 
      detailed = TRUE,
      paired = TRUE
    ) %>%
    add_xy_position(x = "group") %>% 
    adjust_pvalue() %>% 
    add_significance("p.adj")
  
  return(stat.test)
}

paired_adjusted_t_test <- function(dt_melt, id_colnames, x_aes, output_path) {
  # Pairwise t-test between groups
  stat.test <- dt_melt %>%
    group_by(Dimensions) %>%
    # Pairwise t-test between groups
    # If you want to assume the equality of variances (Student t-test), specify the option var.equal = TRUE:
    # https://www.datanovia.com/en/lessons/t-test-in-r/
    t_test(
      values ~ group, 
      # var.equal = TRUE, 
      detailed = TRUE,
      paired = TRUE
    ) %>%
    adjust_pvalue(method = 'bonferroni') %>% 
    add_significance("p.adj")
  
  stat.df <- data.frame(stat.test)
  # Write the dataframe to a CSV file
  write.csv(stat.df, output_path, row.names = FALSE)
  return(as.data.frame(stat.df))
}

paired_adjusted_rstatix_wilcox_test <- function(dt_melt, id_colnames, x_aes, output_path) {
  # Pairwise t-test between groups
  stat.test <- dt_melt %>%
    group_by(Dimensions) %>%
    # Pairwise t-test between groups
    # If you want to assume the equality of variances (Student t-test), specify the option var.equal = TRUE:
    # https://www.datanovia.com/en/lessons/t-test-in-r/
    rstatix::wilcox_test(
      values ~ group, 
      # var.equal = TRUE, 
      detailed = TRUE,
      paired = TRUE
    ) %>%
    adjust_pvalue(method = 'bonferroni') %>% 
    add_significance("p.adj")
  stat.df <- data.frame(stat.test)
  # Write the dataframe to a CSV file
  write.csv(stat.df, output_path, row.names = FALSE)
  return(as.data.frame(stat.df))
}

rstatix_t_test <- function(df) {
  stat.test <- df %>%
    # group_by(group) %>%
    t_test(
      exp ~ group, 
      # var.equal = TRUE, 
      detailed = TRUE) %>%
    add_xy_position(x = "group")
  return(stat.test)
}

rstatix_wilcox_test <- function(df) {
  stat.test <- df %>%
    # group_by(group) %>%
    rstatix::wilcox_test(
      exp ~ group, 
      # var.equal = TRUE, 
      detailed = TRUE) %>%
    add_xy_position(x = "group")
  return(stat.test)
}
# 
# adjusted_t_test <- function(dt_melt, id_colnames, x_aes, output_path) {
#   # Pairwise t-test between groups
#   stat.test <- dt_melt %>%
#     group_by(Dimensions) %>%
#     # Pairwise t-test between groups
#     # If you want to assume the equality of variances (Student t-test), specify the option var.equal = TRUE:
#     # https://www.datanovia.com/en/lessons/t-test-in-r/
#     t_test(
#       values ~ group, 
#       # var.equal = TRUE, 
#       detailed = TRUE) %>%
#     adjust_pvalue(method = 'bonferroni')
#   stat.df <- data.frame(stat.test)
#   # Write the dataframe to a CSV file
#   write.csv(stat.df, output_path, row.names = FALSE)
#   return(as.data.frame(stat.df))
# }

adjusted_t_test <- function(dt_melt, id_colnames, x_aes, output_path, adjust_method = "bonferroni") {
  # Pairwise t-test between groups
  stat.test <- dt_melt %>%
    group_by(Dimensions) %>%
    # Pairwise t-test between groups
    # If you want to assume the equality of variances (Student t-test), specify the option var.equal = TRUE:
    # https://www.datanovia.com/en/lessons/t-test-in-r/
    t_test(
      values ~ group, 
      # var.equal = TRUE, 
      detailed = TRUE) %>%
    adjust_pvalue(method = adjust_method) %>% 
    add_significance("p.adj")
  stat.df <- data.frame(stat.test)
  # Write the dataframe to a CSV file
  write.csv(stat.df, output_path, row.names = FALSE)
  return(as.data.frame(stat.df))
}

adjusted_wilcox_test <- function(dt_melt, id_colnames, x_aes, output_path, adjust_method = "bonferroni") {
  # Pairwise t-test between groups
  stat.test <- dt_melt %>%
    group_by(Dimensions) %>%
    # Pairwise t-test between groups
    # If you want to assume the equality of variances (Student t-test), specify the option var.equal = TRUE:
    # https://www.datanovia.com/en/lessons/t-test-in-r/
    rstatix::wilcox_test(
      values ~ group, 
      # var.equal = TRUE, 
      detailed = TRUE) %>%
    adjust_pvalue(method = adjust_method) %>% 
    add_significance("p.adj")
  stat.df <- data.frame(stat.test)
  # Write the dataframe to a CSV file
  write.csv(stat.df, output_path, row.names = FALSE)
  return(as.data.frame(stat.df))
}

adjusted_paired_boxplot_nojitter <- function(input_df, id_colnames, x_aes, plot_dir, plot_suffix, comparisons_list, plot_width, ncol){
  make_dir(plot_dir)
  id_colname <- setdiff(id_colnames, c(x_aes))
  print(id_colname)
  
  df <- input_df
  df[[id_colname]] <- gsub("-Ca.qptiff|-ca.qptiff|-m.qptiff|-n.qptiff|-M.qptiff|-N.qptiff", "", df[[id_colname]])
  df[[id_colname]] <- gsub("TLS\\s*\\d+\\s*", "", df[[id_colname]])
  
  names(df)[names(df) == x_aes] <- "group"
  # Replace "banana" with "grape"
  id_colnames[id_colnames == x_aes] <- "group"
  plot_path <- paste0(plot_dir, plot_suffix)
  dt_melt <- melt(data.table(df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
  # nomelt_boxplot_nojitter(df, dt_melt, id_colnames, "group", plot_path, comparisons_list, plot_width, ncol)
  stat_df <- paired_adjusted_rstatix_wilcox_test(dt_melt, id_colnames, x_aes, paste0(plot_dir, 'adjusted_p.csv'))
  
  significant_stat_df <- subset(stat_df, p < 0.05)
  significant_stat_df$Dimensions <- as.character(significant_stat_df$Dimensions)
  significant_dimensions <- as.vector(unique(significant_stat_df$Dimensions))
  
  insignificant_stat_df <- subset(stat_df, p >= 0.05)
  insignificant_stat_df$Dimensions <- as.character(insignificant_stat_df$Dimensions)
  insignificant_dimensions <- as.vector(unique(insignificant_stat_df$Dimensions))
  
  single_paired_boxplot_rstatix(
    filter_df(c(id_colnames, significant_dimensions), df),
    id_colnames, "group", paste0(plot_dir, 'significant/'), plot_suffix, comparisons_list
  )
  single_paired_boxplot_rstatix(
    filter_df(c(id_colnames, insignificant_dimensions), df),
    id_colnames, "group", paste0(plot_dir, 'insignificant/'), plot_suffix, comparisons_list
  )
}

single_paired_boxplot_rstatix <- function(df, id_colnames, x_aes, plot_dir, plot_suffix, comparisons_list){
  make_dir(plot_dir)
  suppressWarnings({
    exp_colnames <- setdiff(colnames(df), id_colnames)
    id_colname <- setdiff(id_colnames, c(x_aes))
    for (column in exp_colnames) {
      print(column[1])
      exp_df <- as.data.frame(filter_df(c(id_colnames, column), df))
      names(exp_df)[names(exp_df) == column] <- "exp"
      for (sites in comparisons_list) {
        
        data_df <- as.data.frame(exp_df %>% 
                                   filter(group %in% sites))
        concatenated_sites <- paste(sites, collapse = "_")
        
        wilcox.stat = paired_rstatix_wilcox_test(
          # melt(data.table(exp_df), variable.name="group", value.name = "exp")
          data_df
        )
        
        p <- ggpaired(data = data_df,
                      x="group",
                      y="exp",
                      color="group",
                      ylab = "values",
                      xlab = "group",
                      line.color = "gray", 
                      line.size = 0.4,
                      id = id_colname, 
                      short.panel.labs = FALSE,
                      order = sites
        ) +
          stat_pvalue_manual(wilcox.stat, label = "p", tip.length = 0)
        #  labs(subtitle = get_test_label(t.stat, detailed = TRUE))
        column <- gsub(':', '', column)
        column <- gsub('%', '', column)
        ggsave(
          paste0(plot_dir, column, '.', concatenated_sites, '.', plot_suffix), 
          # width = plot_width, 
          # height = plot_width / ncol * 1.5 * plot_rownum,  # as.integer(float_value)
          dpi = 600,
          limitsize = FALSE
        )
      }
    }
  })
}

adjusted_paired_wilcox_test <- function(df, id_colnames, x_aes, output_path, comparisons_list) {
  for (comparison in comparisons_list) {
    filtered_df <- 
      dt_melt <- melt(data.table(filtered_df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
  }
  # Pairwise t-test between groups
  stat.test <- dt_melt %>%
    group_by(Dimensions) %>%
    do(wilcox = wilcox.test(values ~ group, data = ., paired = TRUE)) %>%
    summarise(id1, Wilcox = wilcox$p.value)
  stat.df <- data.frame(stat.test)
  # Write the dataframe to a CSV file
  write.csv(stat.df, output_path, row.names = FALSE)
}

paired_wilcox_test <- function(df, groups, id_cols) {
  suppressWarnings({
    all_colnames <- colnames(df)
    values_colnames <- all_colnames[!all_colnames %in% id_cols]
    pair1 = groups[1]
    pair2 = groups[2]
    dimensions = c()
    pair_info = c()
    wilcox_p_value = c()
    for (colname in values_colnames) {
      pair1_keyword <- paste0('_',pair1,'_')
      pair2_keyword <- paste0('_',pair2,'_')
      if (grepl(pair1_keyword, colname)) {
        pair2_colname <- gsub(pair1_keyword, pair2_keyword, colname)
        bare_colname <- gsub(pair1_keyword, '', colname)
        # print(bare_colname)
        stat.test <- wilcox.test(df[[colname]], df[[pair2_colname]], paired = TRUE)
        #summarise(bare_colname, pair1_keyword, pair2_keyword, Wilcox = wilcox$p.value)
        dimensions <- c(dimensions, bare_colname)
        pair_info <- c(pair_info, paste0(pair1,'_', pair2))
        wilcox_p_value <- c(wilcox_p_value, stat.test$p.value)
      }
    }
    # Create DataFrame
    df <- data.frame(dimensions = dimensions, pair_info = pair_info, wilcox_p_value = wilcox_p_value)
    return(df)
  })
}

ttest_unit_test <- function(df, id_colnames, x_aes) {
  dt_melt <- melt(data.table(df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
  
  df_wilcox <- df %>%
    group_by(Dimensions) %>%
    do(wilcox_test = tidy(wilcox.test(value.name ~ x_aes, data = ., paired = TRUE)))
  
  
  df_ttest <- dt_melt %>%
    group_by(Dimensions) %>%
    do(t_test = tidy(t.test(values ~ .data[[x_aes]], data = .)))
  print(df_ttest)
}

ttest_adjusted <- function(df, id_colnames, x_aes, plot_path, comparisons_list, plot_width, ncol) {
  suppressWarnings({
    dt_melt <- melt(data.table(df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
    
    data_colnum <- ncol(df) - 1  # Exclude the ID column
    plot_rownum <- data_colnum / ncol
    p <- ggplot(data = dt_melt,
                # aes_string(my_aes)
                aes(
                  x=.data[[x_aes]], 
                  y=values, 
                  color=.data[[x_aes]]
                )
    ) + geom_boxplot() + 
      theme_classic2() +
      theme(legend.position = "none") +  # no color legend 
      #scale_fill_jco() + 
      facet_wrap(~ Dimensions, scales="free", ncol = ncol) +
      # facet_grid(~ Dimensions, scales="free", , space='free') +
      stat_compare_means(aes(method = "t.test", label = "p = {p.format}"), comparisons = comparisons_list)
    pg <- ggplot_build(p)
    print(pg)
  })
}

complex_boxplot_nojitter <- function(
  df, 
  id_colnames, 
  x_aes,
  color_aes,
  shape_aes,
  plot_path, 
  comparisons_list, 
  plot_width, 
  ncol
){
  suppressWarnings({
    dt_melt <- melt(data.table(df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
    data_colnum <- ncol(df) - length(id_colnames)  # Exclude the ID column
    plot_rownum <- data_colnum / ncol
    p <- ggplot(data = dt_melt,
                aes(
                  x=.data[[x_aes]], 
                  y=values, 
                  color=.data[[color_aes]],
                  shape=.data[[shape_aes]]
                )
    ) + geom_boxplot(
      position=position_dodge(0.8)
    ) + 
      theme_classic2() + 
      #scale_fill_jco() + 
      facet_wrap(~ Dimensions, scales = "free", ncol = ncol) +
      stat_compare_means(
        aes(
          method = "t.test", 
          #label = "p = {p.format}",
          label = "p.format",
          paired = TRUE
        )
        , comparisons = comparisons_list
      )
    ggsave(
      plot_path,
      width = plot_width, 
      height = plot_width / ncol * 1.5 * plot_rownum,  # as.integer(float_value)
      limitsize = FALSE)
  })
}

complex_boxplot <- function(
  df, 
  id_colnames, 
  x_aes,
  color_aes,
  shape_aes,
  plot_path, 
  comparisons_list, 
  plot_width, 
  plot_height
){
  dt_melt <- melt(data.table(df), id.vars=id_colnames, variable.name="Dimensions", value.name = "values")
  p <- ggplot(data = dt_melt,
              aes(
                x=.data[[x_aes]], 
                y=values, 
                color=.data[[color_aes]],
                shape=.data[[shape_aes]]
              )
  ) + geom_boxplot(
    position=position_dodge(0.8)
  ) + 
    theme_classic2() + 
    #scale_fill_jco() + 
    facet_wrap(~ Dimensions, ncol = 6, scales = "free") +
    stat_compare_means(
      aes(
        method = "t.test", 
        #label = "p = {p.format}",
        label = "p.format",
        paired = TRUE
      )
      , comparisons = comparisons_list
    ) +
    geom_jitter(shape=16, 
                # position=position_jitter(0.2)
                position=position_dodge(0.8)
    )
  ggsave(plot_path, width = plot_width, height = plot_height, limitsize = FALSE)
}

p_barplot <- function(plot_path, df, colname, group_name, group_order){
  df$Class <- ifelse(df[[colname]] >= median(df[[colname]]), "High", "Low")
  df$Class <- factor(df$Class, levels = c("High", "Low"))
  pdf(file = plot_path, width = 6, height = 6, onefile = F)
  pvalue_bar(data = df, Group="Class", Response=group_name, 
             Group_order=c("High", "Low"), Response_order=group_order, 
             col = c("#2166AC","#66C2A5","#FDBF6F","#B3B3B3"), 
             test_Group = list(c("High", "Low")))
  dev.off()
}

merge_group <- function(group_info_df, df, group_identifier_colname, identifier_colname){
  df[, group_identifier_colname] <- df[, identifier_colname]
  # Merge the dataframes based on the 'ID' column
  merged_df <- left_join(df, group_info_df, by = group_identifier_colname)
  # Delete column from the dataframe
  if (group_identifier_colname != identifier_colname){
    merged_df <- merged_df %>% dplyr::select(-c(group_identifier_colname))
  }
  return(merged_df)
}

my_aes_string <- function(key_colname){
  string <- paste0('x=', key_colname, ', y=values, color=Dimensions')
  return(string)
}

# heatmap function
simple_hp <- function(df, hp_anno, cluster_column){
  hp <- Heatmap(
    t(as.matrix(scale_matrix(numeric_only_df(df)))),
    show_column_names = TRUE, 
    show_row_names = TRUE,
    #row_names_side = "left",
    #col = rainbow(4),
    column_names_rot = -90,
    row_names_gp =  gpar(fontsize = 8),
    column_names_gp =  gpar(fontsize = 8),
    cluster_rows = T,
    clustering_method_rows = "ward.D2",
    name = "Scale(exp)",
    cluster_column_slices = T,
    cluster_columns = cluster_column,
    top_annotation = hp_anno
  )
  print(hp)
  return(hp)
}

pick_color <- function(data_df, colnam, brewer_set){
  selected_column <- data_df[[colnam]]
  mycolors = c(brewer.pal(length(unique(selected_column)), brewer_set)) # Set2
  anno_colors = mycolors[1:length(unique(selected_column))]
  names(anno_colors) = unique(selected_column)
  return(anno_colors)
}

sort_df <- function(df, key_colname){
  # Sort the dataframe by the specified column in ascending order
  sorted_df <- df[order(df[[key_colname]]), ]
  return(sorted_df)
}