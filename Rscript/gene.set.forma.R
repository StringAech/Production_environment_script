library(clusterProfiler)
library(dplyr)
library(openxlsx)
library(magrittr)
library(stringr)
format_gene_set <- function(data,patter="IFN",type=c("list","data.frame")) {
  #### data: the MSigDB geneset list or other database geneset list 
  #### patter: the character you want to retrieve
  #### type: the results you want to get
  message(paste0("format the gene set frome ",patter))
  message(paste0("The results type is ",type))
  IFN_names <- grep(pattern = patter,x = names(c2_sets),value = T)
  IFN_pathway <- c2_sets[IFN_names]
  IFN_pathway_data.frame <- lapply(names(IFN_pathway), function(.x){
    data <- data.frame(
      pathway = .x,
      gene = IFN_pathway[[.x]]
    )
    return(data)
  })
  names(IFN_pathway_data.frame) <- names(IFN_pathway)
  if(type == "list") { return(IFN_pathway_data.frame) }
  if(type == "data.frame") {return(Reduce(rbind,IFN_pathway_data.frame) %>% as.data.frame())}
}
c2_sets <- readRDS("../MsigDB/Hs.c2.all.v7.1.entrez.rds")
IFN_pathway <- format_gene_set(c2_sets,patter = "IFN",type = "data.frame")
CD8_pathway <- format_gene_set(c2_sets,patter = "CD8",type = "data.frame")
