#### add confidence_interval ####
setMethod("PCA_plot", "EPARS_data", function(object, group = c("SegmentLabel"),
                                             CachePath = CachePath, template = template,
                                             color_palette = brewer.pal(name="Accent", n = 8),
                                             outlier_label = F, interactive = F,
                                             text =paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel, '<br>Scan_ID:', Scan_ID),
                                             space = space, FigWidth = 8, FigHeight = 7, prefix = "PCA"){
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggsci))
  df_pca <- prcomp(as.data.frame(t(object@exp)))
  df_pcs <- cbind(df_pca$x, object@anno)
  percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
  percentage <- paste(colnames(df_pcs),"(",paste(as.character(percentage),"%",")",sep=""))
  
  if (interactive == F){
    # when interactive == F which is default, opt static dr plot
    for (i in group) {
      
      p <- ggplot(df_pcs, aes(x = PC1, y = PC2, color = as.character(df_pcs[[i]]))) +
        geom_point(size = 3, alpha = 0.8) + theme_bw() +
        xlab(percentage[1]) + scale_color_manual(values = color_palette) +
        labs(x = percentage[1], y = percentage[2], color = i)+
        geom_hline(yintercept = 0, linetype = 3, size = 0.3) +
        geom_vline(xintercept = 0, linetype = 3, size = 0.3) +
        stat_ellipse(level = 0.95, show.legend = F)
      # added section for outlier label
      if(outlier_label == T){
        p <- p + ggrepel::geom_label_repel(aes(label = rownames(df_pcs)), data = df_pcs,
                                           color = "black", max.overlaps = 2,
                                           box.padding = 1, label.size = 0.06)
      }
      cat(sprintf(template, paste0(i)))
      #print(htmltools::tagList(plotly::ggplotly(p)))
      subchunkify(p, fig_width = FigWidth, fig_height = FigHeight)
      cat(space)
      ggsave(paste0(PlotsPath, prefix, "_", i, ".pdf"), plot = p, width = FigWidth, height = FigHeight)
    }
  } else {
    for (i in group) {
      
      p <-paste0( "plot_ly(df_pcs, x = ~PC1, y = ~PC2, color = ~get(eval(i)),
                   colors = color_palette[1:length(unique(sort(df_pcs[[i]])))],
                   width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
                   marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
                   text =~",  rlang::enquo(text)%>% rlang::as_label(),")") %>% parse(text = .) %>% eval
      p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = percentage[1],
                                                                    zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
                        yaxis = list(title = percentage[2], zeroline = FALSE))
      
      # p <- plot_ly(df_pcs, x = ~PC1, y = ~PC2, color = ~get(eval(i)),
      #              colors = color_palette[1:length(unique(sort(df_pcs[[i]])))],
      #              width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
      #              marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
      #              text = ~paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel, '<br>Scan_ID:', Scan_ID))
      # p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = percentage[1],
      #                   zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
      #                   yaxis = list(title = percentage[2], zeroline = FALSE))
      cat(sprintf(template, paste0(i)))
      print(htmltools::tagList(p))
      cat(space)
    }
  }
  cat(substring(gsub("%s ", "", template, perl = TRUE), 2))
  
  object@dr_result$PCA <- as.matrix(df_pcs)
  write.table(df_pcs, file = paste0(CachePath, prefix, ".csv"), sep = ",", row.names = T, col.names = NA)
  return(object)
})
setMethod("UMAP_plot", "EPARS_data", function(object, group = c("SegmentLabel"), need_log = TRUE,
                                              CachePath = CachePath, template = template,
                                              color_palette = brewer.pal(name="Accent", n = 8),
                                              outlier_label = F, n_neighbors = 15, interactive = F,
                                              text =paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel, '<br>Scan_ID:', Scan_ID),
                                              space = space, FigWidth = 8, FigHeight = 7, prefix = "UMAP"){
  suppressPackageStartupMessages(library(umap))
  if (need_log == TRUE){umap_result <- umap(t(log2(object@exp)),
                                            n_neighbors = n_neighbors, random_state = 2713)}
  if (need_log == FALSE){umap_result <- umap(t(object@exp),
                                             n_neighbors = n_neighbors, random_state = 2713)}
  
  colnames(umap_result$layout) <- c("UMAP_Dim_1","UMAP_Dim_2")
  df_umap <- cbind(umap_result$layout, object@anno)
  
  if (interactive == F){
    # when interactive == F which is default, opt static dr plot
    for (i in group){
      p <- ggplot(df_umap,aes(x=UMAP_Dim_1,y=UMAP_Dim_2,color= as.character(df_umap[[i]])))+
        geom_point(size = 3, alpha = 0.8) + theme_bw() +
        labs(color = i)+ scale_color_manual(values = color_palette) +
        theme(aspect.ratio = 1)+
        stat_ellipse(level = 0.95, show.legend = F)
      # geom_mark_hull(aes(fill = df_umap[[i]]), expand = unit(3, "mm")) +
      # geom_voronoi_segment()
      # added section for outlier label
      if(outlier_label == T){
        p <- p + ggrepel::geom_label_repel(aes(label = rownames(df_umap)), data = df_umap,
                                           color = "black", max.overlaps = 2,
                                           box.padding = 1, label.size = 0.06)
      }
      cat(sprintf(template, paste0(i)))
      subchunkify(p, fig_width = FigWidth, fig_height = FigHeight)
      cat(space)
      ggsave(paste0(PlotsPath, prefix, "_", i, ".pdf"), plot = p, width = FigWidth, height = FigHeight)
    }
  } else {
    for (i in group) {
      p <-paste0( "plot_ly(df_umap, x = ~UMAP_Dim_1, y = ~UMAP_Dim_2, color = ~get(eval(i)),
                   colors = color_palette[1:length(unique(sort(df_umap[[i]])))],
                   width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
                   marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
                   text =~",  rlang::enquo(text)%>% rlang::as_label(),")") %>% parse(text = .) %>% eval
      p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = 'UMAP_Dim_1',
                                                                    zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
                        yaxis = list(title = 'UMAP_Dim_2', zeroline = FALSE))
      
      # p <- plot_ly(df_umap, x = ~UMAP_Dim_1, y = ~UMAP_Dim_2, color = ~get(eval(i)),
      #              colors = color_palette[1:length(unique(sort(df_umap[[i]])))],
      #              width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
      #              marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
      #              text = ~paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel,
      #                            '<br>Scan_ID:', Scan_ID))
      # p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = 'UMAP_Dim_1',
      #                   zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
      #                   yaxis = list(title = 'UMAP_Dim_2', zeroline = FALSE))
      cat(sprintf(template, paste0(i)))
      print(htmltools::tagList(p))
      cat(space)
    }
  }
  
  cat(substring(gsub("%s ", "", template, perl = TRUE), 2))
  
  object@dr_result$UMAP <- as.matrix(df_umap)
  write.table(df_umap, file = paste0(CachePath, prefix, ".csv"), sep = ",", row.names = T, col.names = NA)
  return(object)
})
setMethod("tSNE_plot", "EPARS_data", function(object, group = c("SegmentLabel"), need_log = TRUE,
                                              CachePath = CachePath, template = template,
                                              color_palette = brewer.pal(name="Accent", n = 8),
                                              outlier_label = F, interactive = F,
                                              text =paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel, '<br>Scan_ID:', Scan_ID),
                                              space = space, FigWidth = 8, FigHeight = 7, prefix = "tSNE"){
  suppressPackageStartupMessages(library(Rtsne))
  tSNE_perplexity <- floor((ncol(object@exp)-1)*.2)
  if (tSNE_perplexity < 2){tSNE_perplexity <- 2}
  set.seed(seed = 2713)
  if (need_log == TRUE){tSNE_result <- Rtsne(t(log2(object@exp)),
                                             perplexity = tSNE_perplexity,
                                             check_duplicates = FALSE)}
  if (need_log == FALSE){tSNE_result <- Rtsne(t(object@exp),
                                              perplexity = tSNE_perplexity,
                                              check_duplicates = FALSE)}
  
  colnames(tSNE_result$Y) <- c("tSNE_Dim_1","tSNE_Dim_2")
  df_tSNE <- cbind(tSNE_result$Y, object@anno)
  
  if (interactive == F){
    # when interactive == F which is default, opt static dr plot
    for (i in group){
      p <- ggplot(df_tSNE,aes(x=tSNE_Dim_1,y=tSNE_Dim_2,color= as.character(df_tSNE[[i]])))+
        geom_point(size = 3, alpha = 0.8) + theme_bw() +
        labs(color = i)+ scale_color_manual(values = color_palette) +
        theme(aspect.ratio = 1)+
        stat_ellipse(level = 0.95, show.legend = F)
      # added section for outlier label
      if(outlier_label == T){
        p <- p + ggrepel::geom_label_repel(aes(label = rownames(df_tSNE)), data = df_tSNE,
                                           color = "black", max.overlaps = 2,
                                           box.padding = 1, label.size = 0.06)
      }
      cat(sprintf(template, paste0(i)))
      subchunkify(p, fig_width = FigWidth, fig_height = FigHeight)
      cat(space)
      ggsave(paste0(PlotsPath, prefix, "_", i, ".pdf"), plot = p, width = FigWidth, height = FigHeight)
    }
  } else {
    for (i in group) {
      p <-paste0( "plot_ly(df_tSNE, x = ~tSNE_Dim_1, y = ~tSNE_Dim_2, color = ~get(eval(i)),
                   colors = color_palette[1:length(unique(sort(df_tSNE[[i]])))],
                   width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
                   marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
                   text =~",  rlang::enquo(text)%>% rlang::as_label(),")") %>% parse(text = .) %>% eval
      p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = 'tSNE_Dim_1',
                                                                    zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
                        yaxis = list(title = 'tSNE_Dim_2', zeroline = FALSE))       
      
      # p <- plot_ly(df_tSNE, x = ~tSNE_Dim_1, y = ~tSNE_Dim_2, color = ~get(eval(i)),
      #              colors = color_palette[1:length(unique(sort(df_tSNE[[i]])))],
      #              width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
      #              marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
      #              text = ~paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel,
      #                            '<br>Scan_ID:', Scan_ID))
      # p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = 'tSNE_Dim_1',
      #                   zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
      #                   yaxis = list(title = 'tSNE_Dim_2', zeroline = FALSE))
      cat(sprintf(template, paste0(i)))
      print(htmltools::tagList(p))
      cat(space)
    }
  }
  
  cat(substring(gsub("%s ", "", template, perl = TRUE), 2))
  
  object@dr_result$tSNE <- as.matrix(df_tSNE)
  write.table(df_tSNE, file = paste0(CachePath, prefix, ".csv"), sep = ",", row.names = T, col.names = NA)
  return(object)
})

#### del the "Other" sample ####
    setMethod("PCA_plot", "EPARS_data", function(object, group = c("SegmentLabel"),
                                                 CachePath = CachePath, template = template,
                                                 color_palette = brewer.pal(name="Accent", n = 8),
                                                 outlier_label = F, interactive = F,
                                                 text =paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel, '<br>Scan_ID:', Scan_ID),
                                                 space = space, FigWidth = 8, FigHeight = 7, prefix = "PCA"){
      suppressPackageStartupMessages(library(ggplot2))
      suppressPackageStartupMessages(library(ggsci))
      df_pca <- prcomp(as.data.frame(t(object@exp)))
      df_pcs <- cbind(df_pca$x, object@anno)
      percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
      percentage<- paste(colnames(df_pcs),"(",paste(as.character(percentage),"%",")",sep=""))
      
      if (interactive == F){
        # when interactive == F which is default, opt static dr plot
        for (i in group) {
          
          p <- ggplot(df_pcs, aes(x = PC1, y = PC2, color = as.character(df_pcs[[i]]))) +
            geom_point(size = 3, alpha = 0.8) + theme_bw() +
            xlab(percentage[1]) + scale_color_manual(values = color_palette) +
            labs(x = percentage[1], y = percentage[2], color = i)+
            geom_hline(yintercept = 0, linetype = 3, size = 0.3) +
            geom_vline(xintercept = 0, linetype = 3, size = 0.3) # +
          # stat_ellipse(level = 0.95, show.legend = F)
          # added section for outlier label
          if(outlier_label == T){
            p <- p + ggrepel::geom_label_repel(aes(label = rownames(df_pcs)), data = df_pcs,
                                               color = "black", max.overlaps = 2,
                                               box.padding = 1, label.size = 0.06)
          }
          cat(sprintf(template, paste0(i)))
          #print(htmltools::tagList(plotly::ggplotly(p)))
          subchunkify(p, fig_width = FigWidth, fig_height = FigHeight)
          cat(space)
          ggsave(paste0(PlotsPath, prefix, "_", i, ".pdf"), plot = p, width = FigWidth, height = FigHeight)
        }
      } else {
        for (i in group) {
          
          p <-paste0( "plot_ly(df_pcs, x = ~PC1, y = ~PC2, color = ~get(eval(i)),
                   colors = color_palette[1:length(unique(sort(df_pcs[[i]])))],
                   width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
                   marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
                   text =~",  rlang::enquo(text)%>% rlang::as_label(),")") %>% parse(text = .) %>% eval
          p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = percentage[1],
                                                                        zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
                            yaxis = list(title = percentage[2], zeroline = FALSE))
          
          # p <- plot_ly(df_pcs, x = ~PC1, y = ~PC2, color = ~get(eval(i)),
          #              colors = color_palette[1:length(unique(sort(df_pcs[[i]])))],
          #              width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
          #              marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
          #              text = ~paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel, '<br>Scan_ID:', Scan_ID))
          # p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = percentage[1],
          #                   zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
          #                   yaxis = list(title = percentage[2], zeroline = FALSE))
          cat(sprintf(template, paste0(i)))
          print(htmltools::tagList(p))
          cat(space)
        }
      }
      cat(substring(gsub("%s ", "", template, perl = TRUE), 2))
      
      object@dr_result$PCA <- as.matrix(df_pcs)
      write.table(df_pcs, file = paste0(CachePath, prefix, "_", group, ".csv"), sep = ",", row.names = T, col.names = NA)
      return(object)
    })
    
    # PCA plot
    template <- "#### %s {.unlisted .unnumbered}
"
    for (i in exp_group_object@heatmap_groups){
      exp_group_object_single = exp_group_object
      exp_group_object_single@heatmap_groups = c(i)
      exp_group_object_single@anno <- exp_group_object_single@anno %>% 
        dplyr::filter(.[[i]] != "Other")
      exp_group_object_single@exp <- exp_group_object_single@exp[,rownames(exp_group_object_single@anno)]
      exp_group_object_single <- PCA_plot(exp_group_object_single, CachePath = CachePath, 
                                          group = exp_group_object_single@heatmap_groups, 
                                          color_palette = exp_group_object@config$pan_color,
                                          template = template, space = space,
                                          prefix = "PCA")
    }
    
    ### tSNE降维 {.tabset .tabset-fade}
    setMethod("tSNE_plot", "EPARS_data", function(object, group = c("SegmentLabel"), need_log = TRUE,
                                                  CachePath = CachePath, template = template,
                                                  color_palette = brewer.pal(name="Accent", n = 8),
                                                  outlier_label = F, interactive = F,
                                                  text =paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel, '<br>Scan_ID:', Scan_ID),
                                                  space = space, FigWidth = 8, FigHeight = 7, prefix = "tSNE"){
      suppressPackageStartupMessages(library(Rtsne))
      tSNE_perplexity <- floor((ncol(object@exp)-1)*.2)
      if (tSNE_perplexity < 2){tSNE_perplexity <- 1}
      set.seed(seed = 2713)
      if (need_log == TRUE){tSNE_result <- Rtsne(t(log2(object@exp)),
                                                 perplexity = tSNE_perplexity,
                                                 check_duplicates = FALSE)}
      if (need_log == FALSE){tSNE_result <- Rtsne(t(object@exp),
                                                  perplexity = tSNE_perplexity,
                                                  check_duplicates = FALSE)}
      
      colnames(tSNE_result$Y) <- c("tSNE_Dim_1","tSNE_Dim_2")
      df_tSNE <- cbind(tSNE_result$Y, object@anno)
      
      if (interactive == F){
        # when interactive == F which is default, opt static dr plot
        for (i in group){
          p <- ggplot(df_tSNE,aes(x=tSNE_Dim_1,y=tSNE_Dim_2,color= as.character(df_tSNE[[i]])))+
            geom_point(size = 3, alpha = 0.8) + theme_bw() +
            labs(color = i)+ scale_color_manual(values = color_palette) +
            theme(aspect.ratio = 1)
          # added section for outlier label
          if(outlier_label == T){
            p <- p + ggrepel::geom_label_repel(aes(label = rownames(df_tSNE)), data = df_tSNE,
                                               color = "black", max.overlaps = 2,
                                               box.padding = 1, label.size = 0.06)
          }
          cat(sprintf(template, paste0(i)))
          subchunkify(p, fig_width = FigWidth, fig_height = FigHeight)
          cat(space)
          ggsave(paste0(PlotsPath, prefix, "_", i, ".pdf"), plot = p, width = FigWidth, height = FigHeight)
        }
      } else {
        for (i in group) {
          p <-paste0( "plot_ly(df_tSNE, x = ~tSNE_Dim_1, y = ~tSNE_Dim_2, color = ~get(eval(i)),
                   colors = color_palette[1:length(unique(sort(df_tSNE[[i]])))],
                   width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
                   marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
                   text =~",  rlang::enquo(text)%>% rlang::as_label(),")") %>% parse(text = .) %>% eval
          p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = 'tSNE_Dim_1',
                                                                        zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
                            yaxis = list(title = 'tSNE_Dim_2', zeroline = FALSE))       
          
          # p <- plot_ly(df_tSNE, x = ~tSNE_Dim_1, y = ~tSNE_Dim_2, color = ~get(eval(i)),
          #              colors = color_palette[1:length(unique(sort(df_tSNE[[i]])))],
          #              width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
          #              marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
          #              text = ~paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel,
          #                            '<br>Scan_ID:', Scan_ID))
          # p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = 'tSNE_Dim_1',
          #                   zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
          #                   yaxis = list(title = 'tSNE_Dim_2', zeroline = FALSE))
          cat(sprintf(template, paste0(i)))
          print(htmltools::tagList(p))
          cat(space)
        }
      }
      
      cat(substring(gsub("%s ", "", template, perl = TRUE), 2))
      
      object@dr_result$tSNE <- as.matrix(df_tSNE)
      write.table(df_tSNE, file = paste0(CachePath, prefix, "_", group, ".csv"), sep = ",", row.names = T, col.names = NA)
      return(object)
    })
    
    # UMAP plot
    template <- "#### %s {.unlisted .unnumbered}
"
    for (i in exp_group_object@heatmap_groups){
      exp_group_object_single = exp_group_object
      exp_group_object_single@heatmap_groups = c(i)
      exp_group_object_single@anno <- exp_group_object_single@anno %>% 
        dplyr::filter(.[[i]] != "Other")
      exp_group_object_single@exp <- exp_group_object_single@exp[,rownames(exp_group_object_single@anno)]
      exp_group_object_single <- tSNE_plot(exp_group_object_single, CachePath = CachePath,
                                           group = exp_group_object_single@heatmap_groups, 
                                           color_palette = exp_group_object_single@config$pan_color,
                                           template = template, space = space,
                                           prefix = "tSNE", need_log = FALSE)
    }

    ### UMAP降维 {.tabset .tabset-fade}
    setMethod("UMAP_plot", "EPARS_data", function(object, group = c("SegmentLabel"), need_log = TRUE,
                                                  CachePath = CachePath, template = template,
                                                  color_palette = brewer.pal(name="Accent", n = 8),
                                                  outlier_label = F, n_neighbors = 15, interactive = F,
                                                  text =paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel, '<br>Scan_ID:', Scan_ID),
                                                  space = space, FigWidth = 8, FigHeight = 7, prefix = "UMAP"){
      suppressPackageStartupMessages(library(umap))
      if (need_log == TRUE){umap_result <- umap(t(log2(object@exp)),
                                                n_neighbors = n_neighbors, random_state = 2713)}
      if (need_log == FALSE){umap_result <- umap(t(object@exp),
                                                 n_neighbors = n_neighbors, random_state = 2713)}
      
      colnames(umap_result$layout) <- c("UMAP_Dim_1","UMAP_Dim_2")
      df_umap <- cbind(umap_result$layout, object@anno)
      
      if (interactive == F){
        # when interactive == F which is default, opt static dr plot
        for (i in group){
          p <- ggplot(df_umap,aes(x=UMAP_Dim_1,y=UMAP_Dim_2,color= as.character(df_umap[[i]])))+
            geom_point(size = 3, alpha = 0.8) + theme_bw() +
            labs(color = i)+ scale_color_manual(values = color_palette) +
            theme(aspect.ratio = 1)
          # geom_mark_hull(aes(fill = df_umap[[i]]), expand = unit(3, "mm")) +
          # geom_voronoi_segment()
          # added section for outlier label
          if(outlier_label == T){
            p <- p + ggrepel::geom_label_repel(aes(label = rownames(df_umap)), data = df_umap,
                                               color = "black", max.overlaps = 2,
                                               box.padding = 1, label.size = 0.06)
          }
          cat(sprintf(template, paste0(i)))
          subchunkify(p, fig_width = FigWidth, fig_height = FigHeight)
          cat(space)
          ggsave(paste0(PlotsPath, prefix, "_", i, ".pdf"), plot = p, width = FigWidth, height = FigHeight)
        }
      } else {
        for (i in group) {
          p <-paste0( "plot_ly(df_umap, x = ~UMAP_Dim_1, y = ~UMAP_Dim_2, color = ~get(eval(i)),
                   colors = color_palette[1:length(unique(sort(df_umap[[i]])))],
                   width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
                   marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
                   text =~",  rlang::enquo(text)%>% rlang::as_label(),")") %>% parse(text = .) %>% eval
          p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = 'UMAP_Dim_1',
                                                                        zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
                            yaxis = list(title = 'UMAP_Dim_2', zeroline = FALSE))
          
          # p <- plot_ly(df_umap, x = ~UMAP_Dim_1, y = ~UMAP_Dim_2, color = ~get(eval(i)),
          #              colors = color_palette[1:length(unique(sort(df_umap[[i]])))],
          #              width = 800, height = 650, alpha = 0.8, type = 'scatter', mode = 'markers',
          #              marker = list(symbol = 'circle', sizemode = 'diameter', size=12),
          #              text = ~paste('AOI:', Unique_Sample_ID, '<br>SegmentLabel:', SegmentLabel,
          #                            '<br>Scan_ID:', Scan_ID))
          # p <- p %>% layout(legend = list(x = 1, y = 0.5), xaxis = list(title = 'UMAP_Dim_1',
          #                   zeroline = FALSE, scaleanchor = "y", scaleratio = 1),
          #                   yaxis = list(title = 'UMAP_Dim_2', zeroline = FALSE))
          cat(sprintf(template, paste0(i)))
          print(htmltools::tagList(p))
          cat(space)
        }
      }
      
      cat(substring(gsub("%s ", "", template, perl = TRUE), 2))
      
      object@dr_result$UMAP <- as.matrix(df_umap)
      write.table(df_umap, file = paste0(CachePath, prefix, "_", group, ".csv"), sep = ",", row.names = T, col.names = NA)
      return(object)
    })
    
    # UMAP plot
    template <- "#### %s {.unlisted .unnumbered}
"
    for (i in exp_group_object@heatmap_groups){
      exp_group_object_single = exp_group_object
      exp_group_object_single@heatmap_groups = c(i)
      exp_group_object_single@anno <- exp_group_object_single@anno %>% 
        dplyr::filter(.[[i]] != "Other")
      exp_group_object_single@exp <- exp_group_object_single@exp[,rownames(exp_group_object_single@anno)]
      exp_group_object_single <- UMAP_plot(exp_group_object_single, CachePath = CachePath,n_neighbors = 4, 
                                           group = exp_group_object_single@heatmap_groups, 
                                           color_palette = exp_group_object_single@config$pan_color,
                                           template = template, space = space,
                                           prefix = "UMAP", need_log = FALSE)
    }

    