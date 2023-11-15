library(clusterProfiler)
library(dplyr)
library(openxlsx)
library(magrittr)
library(stringr)
format_gene_set <- function(data,patter="IFN",type=c("list","data.frame"),is.bitr = FALSE,org = "hsa") {
  #### data: the MSigDB geneset list or other database geneset list 
  #### patter: the character you want to retrieve
  #### type: the results you want to get
  message(paste0("format the gene set frome ",patter))
  message(paste0("The results type is ",type))
  names <- grep(pattern = patter,x = names(data),value = T)
  pathway <- data[names]
  pathway.data.frame <- lapply(names(pathway), function(.x){
    data <- data.frame(
      pathway = .x,
      gene = pathway[[.x]]
    )
    return(data)
  })
  names(pathway.data.frame) <- names(pathway)
  if(is.bitr){ 
    suppressPackageStartupMessages(library(clusterProfiler))
    switch (org,
            hsa = {
              suppressPackageStartupMessages(library(org.Hs.eg.db))
              pathway.data.frame <- suppressMessages(lapply(names(pathway.data.frame), function(.x){
                eg = bitr(pathway.data.frame[[.x]]$gene, fromType = "ENTREZID", toType = "SYMBOL", 
                          OrgDb = "org.Hs.eg.db") %>% mutate(pathway = .x) %>% dplyr::select(3,2,1)
              }))
              names(pathway.data.frame) <- names(pathway.data.frame)
            },
            mm = {
              suppressPackageStartupMessages(library(org.Mm.eg.db))
              pathway.data.frame <- suppressMessages(lapply(names(pathway.data.frame), function(.x){
                eg = bitr(pathway.data.frame[[.x]]$gene, fromType = "ENTREZID", toType = "SYMBOL", 
                          OrgDb = "org.Mm.eg.db") %>% mutate(pathway = .x) %>% dplyr::select(3,2,1)
              }))
              names(pathway.data.frame) <- names(pathway.data.frame)
            }
    )
  }
  if(type == "list") {return(pathway.data.frame) }
  if(type == "data.frame") {return(Reduce(rbind,pathway.data.frame) %>% as.data.frame())}
}
c2_sets <- readRDS("~/script/MsigDB/Hs.c2.all.v7.1.entrez.rds")
IFN_pathway <- format_gene_set(c2_sets,patter = "IFN",type = "data.frame",is.bitr = T,org = "hsa")
CD8_pathway <- format_gene_set(c2_sets,patter = "CD8",type = "data.frame",is.bitr = T,org = "has")
