## parallel script
parallel_ORA_enrichment <- function(exp_group_object, CachePath, PlotsPath, core.num = 10, 
                                    ...) {
    enrichKEGGfunction <- function(gene, org = "hsa", pvalueCutoff = 1,
                                   qvalueCutoff = 1) {
        rescale.AsIs <- function(x, ...) {
            dropAsis <- function(x) {
                cls <- class(x)
                structure(x, class = setdiff(cls, "AsIs"))
            }
            scales:::rescale(dropAsis(x), ...)
        }
        if (org == "mm") {
            org <- "mmu"
        }
        library(clusterProfiler)
        library(ggplot2)
        gene <- tryCatch(suppressMessages(tranSymbol(gene, org)), error = function(e) {
        })
        if (class(gene) == "try-error") {
            warning("no gene tranSymbol")
        }
        enr <- try(suppressMessages(enrichKEGG(gene = gene$ENTREZID, organism = org,
                                               pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)))
        if (org == "hsa") {
            enr <- tryCatch(suppressMessages(setReadable(enr, OrgDb = org.Hs.eg.db,
                                                         keyType = "ENTREZID")), error = function(e) {
                                                         })
        }
        if (org == "mm") {
            enr <- tryCatch(suppressMessages(setReadable(enr, OrgDb = org.Mm.eg.db,
                                                         keyType = "ENTREZID")), error = function(e) {
                                                         })
        }
        if (class(enr) == "try-error") {
            enr <- NULL
        }
        enr_res <- tryCatch(as.data.frame(enr), error = function(e) {
            data.frame()
        })
        if (is.null(enr) == TRUE) {
            p <- text_only_ggplot("Not enough differentially expressed genes \n for enrichment analysis")
            p2 <- text_only_ggplot("Not enough differentially expressed genes \n for enrichment analysis")
        } else {
            if (org == "hsa") {
                enr <- setReadable(enr, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
            }
            if (org == "mm") {
                enr <- setReadable(enr, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
            }
            enr_res <- enrichres_id_transform(enr_res, gene)
            enr_plot <- enr
            enr_plot@result$Description <- wrap.labels(enr_plot@result$Description,
                                                       50)
            term_num <- length(enr$Description)
            if (term_num > 20) {
                term_num <- 20
            }
            p <- tryCatch(suppressMessages(clusterProfiler::dotplot(enr_plot,
                                                                    showCategory = term_num, font.size = 7.5) + scale_colour_gradient(high = "royalblue3",
                                                                                                                                      low = "indianred3") + theme(plot.margin = unit(c(1, 2,
                                                                                                                                                                                       1, 2), "lines"), legend.text = element_text(size = 7),
                                                                                                                                                                  legend.title = element_text(size = 7))), error = function(e) {
                                                                                                                                                                      text_only_ggplot("Not enough differentially expressed genes \n mapped to KEGG pathway")
                                                                                                                                                                  })
            p2 <- tryCatch(suppressMessages(clusterProfiler::cnetplot(enr_plot,
                                                                      cex_label_category = 0.7, cex_label_gene = 0.6, showCategory = 5,
                                                                      layout = "kk", colorEdge = TRUE) + theme(legend.position = "none",
                                                                                                               panel.border = element_rect(colour = "grey50", fill = NA,
                                                                                                                                           size = 0.5), plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"))),
                           error = function(e) {
                               text_only_ggplot("Not enough differentially expressed genes \n mapped to KEGG pathway")
                           })
            if (inherits(try(ggplot_build(p), silent = TRUE), "try-error")) {
                p <- text_only_ggplot("Not enough differentially expressed genes \n mapped to KEGG pathway")
            }
            if (inherits(try(ggplot_build(p2), silent = TRUE), "try-error")) {
                p2 <- text_only_ggplot("Not enough differentially expressed genes \n mapped to KEGG pathway")
            }
        }
        width <- 6
        height <- 6
        return_list <- list()
        return_list[[1]] <- enr
        return_list[[2]] <- p
        return_list[[3]] <- c(width, height)
        return_list[[4]] <- enr_res
        return_list[[5]] <- p2
        return(return_list)
    }
    enrich.function <- function(g) {
        tmp_deg <- exp_group_object@deg[[g]]@result_df
        filtered_up <- rownames(tmp_deg %>%
                                    dplyr::filter(tmp_deg$change == "UP"))
        filtered_down <- rownames(tmp_deg %>%
                                      dplyr::filter(tmp_deg$change == "DOWN"))
        GO_UP <- enrichGOfunction(filtered_up, org = exp_group_object@config$org)  #pvalueCutoff  = 1, qvalueCutoff  = 1
        GO_DOWN <- enrichGOfunction(filtered_down, org = exp_group_object@config$org)
        KEGG_UP <- enrichKEGGfunction(filtered_up, org = exp_group_object@config$org)  #pvalueCutoff  = 1, qvalueCutoff  = 1
        KEGG_DOWN <- enrichKEGGfunction(filtered_down, org = exp_group_object@config$org)
        write.table(as.data.frame(GO_UP[[1]]), file = paste0(CachePath,
                                                             gsub("的差异", "", exp_group_object@deg[[g]]@label_string),
                                                             "_GOEnrich_UPGene.csv"), sep = ",", row.names = T, col.names = NA)
        write.table(as.data.frame(GO_DOWN[[1]]), file = paste0(CachePath,
                                                               gsub("的差异", "", exp_group_object@deg[[g]]@label_string),
                                                               "_GOEnrich_DOWNGene.csv"), sep = ",", row.names = T, col.names = NA)
        write.table(as.data.frame(KEGG_UP[[4]]), file = paste0(CachePath,
                                                               gsub("的差异", "", exp_group_object@deg[[g]]@label_string),
                                                               "_KEGGEnrich_UPGene.csv"), sep = ",", row.names = T, col.names = NA)
        write.table(as.data.frame(KEGG_DOWN[[4]]), file = paste0(CachePath,
                                                                 gsub("的差异", "", exp_group_object@deg[[g]]@label_string),
                                                                 "_KEGGEnrich_DOWNGene.csv"), sep = ",", row.names = T, col.names = NA)
        return(list(GO_UP = GO_UP, GO_DOWN = GO_DOWN, KEGG_UP = KEGG_UP,
                    KEGG_DOWN = KEGG_DOWN))
        # exp_group_object@deg[[g]]@deg_enrich$GO_UP <-
        # as.data.frame(GO_UP[[1]])
        # exp_group_object@deg[[g]]@deg_enrich$GO_DOWN <-
        # as.data.frame(GO_DOWN[[1]])
        # exp_group_object@deg[[g]]@deg_enrich$KEGG_UP <-
        # as.data.frame(KEGG_UP[[4]])
        # exp_group_object@deg[[g]]@deg_enrich$KEGG_DOWN <-
        # as.data.frame(KEGG_DOWN[[4]])
    }
    CachePath <- CachePath
    PlotsPath <- PlotsPath
    library(EPARS)
    library(parallel)  #加载并行计算包
    cl <- makeCluster(core.num);  # 初始化cpu集群
    clusterExport(cl, 'exp_group_object');
    # clusterExport(cl, 'CachePath');
    # clusterExport(cl, "space")
    # clusterExport(cl, "PlotsPath");
    # clusterExport(cl, "glm_diff")
    # clusterExport(cl, "rescale.AsIs")
    clusterEvalQ(cl, library(dplyr));  #添加并行计算中用到的包
    clusterEvalQ(cl, library(EPARS)); #添加并行计算中用到的包
    clusterEvalQ(cl, library(RColorBrewer))  #添加并行计算中用到的包
    clusterEvalQ(cl, library(R.utils))  #添加并行计算中用到的包
    clusterEvalQ(cl, library(clusterProfiler))  #添加并行计算中用到的包
    clusterEvalQ(cl, library(ggpubr))  #添加并行计算中用到的包
    enrichment.res <- parLapply(cl, names(exp_group_object@deg), enrich.function)
    parallel::stopCluster(cl)
    names(enrichment.res) <- names(exp_group_object@deg)
    for (i in names(exp_group_object@deg)) {
        exp_group_object@deg[[i]]@deg_enrich <- enrichment.res[[i]]
    }
    return(exp_group_object)
}
#### for example
# exp_group_object <- parallel_ORA_enrichment(exp_group_object, 
#                              core.num = 16,
#                              CachePath = "./test/plot/aaaa/",
#                              PlotsPath = "./test/plot/aaaa/")
#### Extraction of image data
# p1_list <- list()
# p2_list <- list()
# 
# for (g in names(exp_group_object@deg)) {
#     p <- ggarrange(exp_group_object@deg[[g]]@deg_enrich$GO_UP[[2]], exp_group_object@deg[[g]]@deg_enrich$KEGG_UP[[2]], 
#                    exp_group_object@deg[[g]]@deg_enrich$GO_DOWN[[2]], exp_group_object@deg[[g]]@deg_enrich$KEGG_DOWN[[2]], 
#                    labels = c("GO enrich for UP regulated genes", "KEGG enrich for UP regulated genes",
#                               "GO enrich for DOWN regulated genes", "KEGG enrich for DOWN regulated genes"),
#                    ncol = 2, nrow = 2, heights = 13, widths = 12)
#     p1_list[[g]] <- p
# }
# 
# for (g in names(exp_group_object@deg)) {
#     p2 <- ggarrange(exp_group_object@deg[[g]]$GO_UP[[4]], exp_group_object@deg[[g]]$KEGG_UP[[5]], 
#                     exp_group_object@deg[[g]]$GO_DOWN[[4]], exp_group_object@deg[[g]]$KEGG_DOWN[[5]], 
#                     labels = c("GO enrich for UP regulated genes", "KEGG enrich for UP regulated genes",
#                                "GO enrich for DOWN regulated genes", "KEGG enrich for DOWN regulated genes"),
#                     ncol = 2, nrow = 2, heights = 13, widths = 12)
#     p2_list[[g]] <- p2
# }