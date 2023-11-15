library(dplyr)
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
file.list <- list.files("./","*.csv")
### annotation DEG file
purrr::walk(file.list,function(.x){
  DEG <- read.csv(.x)
  DEG_merge <- merge(annotation_geneid,DEG,by = 1)
  write.csv(DEG_merge,file = paste0("./add_geneID/geneID_",.x),row.names = F)
})
### annotation raw/nor expression file
raw_exp <- merge(annotation_geneid,exp_group_object@raw_exp,by.x = "TargetName",by.y = 0)
exp <- merge(annotation_geneid,exp_group_object@exp,by.x = "TargetName",by.y = 0)
write.csv(raw_exp,file = "./add_geneID/raw_exp.csv",row.names = F)
write.csv(exp,file = "./add_geneID/exp.csv",row.names = F)
