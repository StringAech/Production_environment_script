library(dplyr)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
`%!in%` <- Negate(`%in%`)
#### try direct annotation ####
gene = bitr(rownames(a), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db") %>% dplyr::arrange(ENTREZID)
gene <- bitr(
  rownames(exp_group_object@exp),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
) %>% arrange(ENTREZID)
setdiff(rownames(exp_group_object@exp),gene$SYMBOL)
#### use Nanostring annotation file ####
annotation_geneid <- read.csv("~/Reference_Data/WTA_symbol2geneID.csv",sep = ",") %>% 
  dplyr::select(c('TargetName','GeneID'))
file.list <- list.files("./","*.xlsx")
### annotation DEG file
purrr::walk(file.list,function(.x){
  # DEG <- read.csv(.x)
  DEG <- read.xlsx(.x)
  # DEG_merge <- full_join(DEG, annotation_geneid, by = join_by(Gene == TargetName )) %>% 
  #   dplyr::select(1, GeneID, everything())
  DEG_merge <- merge(DEG, annotation_geneid, by = 1, all.x = T, sort = F) %>% 
    dplyr::select(1, GeneID, everything())
    
  # write.csv(DEG_merge,file = paste0("./add_geneID/geneID_",.x),row.names = F)
  write.xlsx(x = DEG_merge, file = paste0("./add_geneID/geneID_", .x), overwrite = T)
})
### annotation raw/nor expression file
raw_exp <- merge(annotation_geneid,exp_group_object@raw_exp,by.x = "TargetName",by.y = 0)
exp <- merge(annotation_geneid,exp_group_object@exp,by.x = "TargetName",by.y = 0)
write.csv(raw_exp,file = "./add_geneID/raw_exp.csv",row.names = F)
write.csv(exp,file = "./add_geneID/exp.csv",row.names = F)
