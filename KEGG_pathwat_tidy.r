library(rjson)
library(jsonlite)
library(tidyverse)
library(dplyr)
KO <- fromJSON(file = "/home/rstudio/user/guohao/Reference_Data/hsa00001.json")#下载并解析JSON文件
KO$name <- NULL
KO <- as.data.frame(KO) %>% 
  unnest(cols = c("children.name","children.children"),names_repair = tidyr_legacy) %>%#重要函数
  unnest(cols = c("children.name","name","children"),names_repair = tidyr_legacy) %>%
  unnest(cols = c("children.name","name","name1","children"),names_repair = tidyr_legacy)
colnames(KO) <- c("L1","L2","L3","KO") 
KO %<>% #整理KEGG ORTHOLOGY
  select(last_col(),everything()) %>%
  separate(col = "KO",sep = ";",into = c("KO","Description")) %>%
  separate(col = "L1",sep = " ",into = c("L1_ID","L1"),extra = "merge") %>%
  filter(!L1_ID %in% c("09180","09190")) %>% #去除BRITE hierarchies和Not Included in Pathway or Brite两大类
  separate(col = "L2",sep = " ",into = c("L2_ID","L2"),extra = "merge") %>%
  separate(col = "L3",sep = " ",into = c("L3_ID","L3"),extra = "merge") %>%
  separate(col = "L3",sep = " \\[PATH:",into = c("L3","PathwayID")) %>%
  mutate(PathwayID=str_remove(PathwayID,pattern = "\\]")) %>%
  separate(col = "KO",sep = " ",into = c("GeneID","GeneSymobol")) %>%
  drop_na()#KEGG ORTHOLOGY等级有缺失的删掉
head(KO)

lipid_pathway <- dplyr::filter(KO,L2 == "Lipid metabolism")
L3_pathway <- split(lipid_pathway,lipid_pathway$L3)
L3_pathway %<>% lapply(., function(.x){
  .x <- list(annotation_matrix = .x,genelist = .x$GeneSymbol)
  return(.x)
})
saveRDS(L3_pathway,file = "~/Reference_Data/KEGG_lipid_pathwat_genelist.rds")