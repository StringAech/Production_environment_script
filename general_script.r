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