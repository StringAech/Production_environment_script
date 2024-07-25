enrichGOfunction = function(gene, mgssize = 10, org = "hsa", ont = "ALL", pvalueCutoff  = 1, qvalueCutoff  = 1){
  gene <-  tryCatch(suppressMessages(tranSymbol(gene, org)), error = function(e){})
  if(class(gene)=="try-error") { warning("no gene tranSymbol")}
  ## 开始富集
  if (org == "hsa") {
    enr <-
      tryCatch(
        clusterProfiler::enrichGO(
          gene$ENTREZID,
          OrgDb = 'org.Hs.eg.db',
          ont = ont,
          pvalueCutoff  = pvalueCutoff,
          qvalueCutoff  = qvalueCutoff,
          readable  = TRUE,
          minGSSize = mgssize
        ),
        error = function(e) {
        }
      )
  }
  if (org == "mm") {
    enr <-
      tryCatch(
        clusterProfiler::enrichGO(
          gene$ENTREZID,
          OrgDb = 'org.Mm.eg.db',
          ont = ont,
          pvalueCutoff  = pvalueCutoff,
          qvalueCutoff  = qvalueCutoff,
          readable  = TRUE,
          minGSSize = mgssize
        ),
        error = function(e) {
        }
      )
  }
  
  return(enr)
}
