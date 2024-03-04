`%!in%` <- Negate(`%in%`)
#### trans2Symbol #### @
trans2symbol <- function(.data, RNASeq_mRNA_org = "hsa", gene = TRUE, transcrit = F, ...) {
  # read in gene id
  if (!("Geneid" %in% colnames(.data))) {
    stop(paste0("you need set the ", deparse(substitute(.data)), " colname 'Geneid'"))
  }
  if (gene) {
    if (RNASeq_mRNA_org == "hsa") {
      gene_map <- read.table("/home/rstudio/user/guohao/script/resource/extdata/annotate/gencode.v37.protein_coding.geneid.genename.txt")
    }
    if (RNASeq_mRNA_org == "mm") {
      gene_map <- read.table("/home/rstudio/user/guohao/script/resource/extdata/annotate/gene_mouse.txt")
    }
    colnames(gene_map) <- c("GeneId", "GeneName")
    gene_map <- gene_map[!duplicated(gene_map[, c("GeneName")]), ]
    merge.file <- merge(gene_map, .data, by.x = "GeneId", by.y = "Geneid")
  }
  
  if (transcrit) {
    # read in trancript id
    anno_file <- "/home/rstudio/user/guohao/script/resource/"
    
    if (RNASeq_mRNA_org == "hsa") {
      transcript_map <- read.csv(file = paste0(anno_file, "gencode.v43.metadata.HGNC"), header = F, sep = "\t")
      colnames(transcript_map) <- c("GeneId", "GeneName", "HGNC_ID")
    }
    
    if (RNASeq_mRNA_org == "mm") {
      transcript_map <- read.csv(file = paste0(anno_file, "gencode.vM32.metadata.MGI"), header = F, sep = "\t")
      colnames(transcript_map) <- c("GeneId", "GeneName", "MGI_ID")
    }
    
    merge.file <- merge(transcript_map, .data, by.x = "GeneId", by.y = "Geneid")
  }
  return(merge.file)
}

#### is.sum_equal ####
# 判断行列和是否为某一个数
# data 需要判读的矩阵
# margin 1 行 2 列
# select_col 选择的列, 要求列是数字
# equal 和
# is.exclude 逻辑值 是否要排除掉符合条件的列
## TRUE 则返回不包含 符合条件的列名 反正则只返回 符合条件的列名

is.sum_equal <- function(data, margin = 2, select_col = "ABC", equal = 0, is.exclude = TRUE){
  res <- apply(data[, select_col], MARGIN = margin, \(.x) isTRUE(all.equal(sum(.x), equal))) 
  if (is.exclude) {return(select_col[!res])}else{return(select_col[res])}
}
## eg
# plot.data <- tidy.anno(immcell_exp, anno, number = i, col = "comment2") %>% 
#   dplyr::select(1:2, is.sum_equal(data = ., margin = 2, select_col = colnames(.)[3:26], equal = 0, is.exclude = T))