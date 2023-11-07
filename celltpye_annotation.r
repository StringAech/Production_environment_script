`%!in%` <- Negate(`%in%`)

Smooth <- c("MYH11","ACTA2")
fibroblast <- c("COL1A1","PDGFRA")
Endothelial <- c("PECAM1","CDH5")
epithelial <- c("KRT5","KRT8","EPCAM")
Tcell <- c("CD3E","CD3G","CD3D")
Bcell <- c("MS4A1","CD79A","CD79B")
mononuclear <- c("Ly79c6")
NKcell <- c("KLRD1","NKG7","NCAM1","FCGR3A")
immue <- c("PTPRC")
Eosinophils <- c("ECP","PRG2","RNASE3","EPX","CCR3","IL5RA","ALOX15")

markers <- list(Smooth=Smooth,Fibroblast=fibroblast,Endothelial=Endothelial,Epithelial=epithelial,
             Tcell=Tcell,Bcell=Bcell,NKcell=NKcell,Eosinophils=Eosinophils)
#### the first annotation ####
p <- DotPlot(Seurat_object, 
             features = c("PTPRC","EPCAM","PECAM1","MME","MYH11","ACTA2","COL1A1","PDGFRA"), 
             group.by = "integrated_snn_res.0.4",)+
  RotatedAxis()+
  scale_x_discrete("")+
  scale_y_discrete("")
ggsave(filename = "./20231105_celltype/the first annotation.pdf",plot = p)


Seurat_object$Celltype_main_level <- "ABC"
attach(Seurat_object@meta.data)
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(4,5,7,14,16,11)),]$Celltype_main_level <- "Immune"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(0,1,2,6,8,10,12,13,15)),]$Celltype_main_level <- "Epithelial"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(17)),]$Celltype_main_level <- "Endothelial"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(9)),]$Celltype_main_level <- "Stromal"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(3)),]$Celltype_main_level <- "Smooth_muscle_cell"
detach(Seurat_object@meta.data)
p <- DotPlot(Seurat_object, features = markers,group.by = "Celltype_main_level",)+
  RotatedAxis()+
  scale_x_discrete("")+
  scale_y_discrete("")
ggsave("./20231105_celltype/the_first_annotation_results.pdf",width = 18,height = 8)
#### Cluster immune ####
immnue_marker <- list(Tcell = Tcell,Bcell = Bcell,NKcell = NKcell)
p <- DotPlot(Seurat_object, features = immnue_marker,group.by = "integrated_snn_res.0.4",)+
  RotatedAxis()+
  scale_x_discrete("")+
  scale_y_discrete("")
ggsave(filename = "./20231105_celltype/immnue cluster annotation.pdf",plot = p,width = 18,height = 8)
attach(Seurat_object@meta.data)
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(7)),]$Celltype_main_level <- "NKcell"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(14,11)),]$Celltype_main_level <- "Bcell"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(4)),]$Celltype_main_level <- "Tcell"
detach(Seurat_object@meta.data)
p <- DotPlot(Seurat_object, features = markers,group.by = "Celltype_main_level",)+
  RotatedAxis()+
  scale_x_discrete("")+
  scale_y_discrete("")
ggsave("./20231105_celltype/immnue cluster annotation result.pdf",width = 18,height = 8)

#### Cluster Other cell type ####
Other_marker <- list(Eosinophils = Eosinophils,Smooth = Smooth,
                     Fibroblast = fibroblast,Endothelial = Endothelial,
                     Epithelial = epithelial)
p <- DotPlot(Seurat_object, features = Other_marker,group.by = "integrated_snn_res.0.4",)+
  RotatedAxis()+
  scale_x_discrete("")+
  scale_y_discrete("")
ggsave(filename = "./20231105_celltype/Other cluster annotation.pdf",plot = p,width = 18,height = 8)
attach(Seurat_object@meta.data)
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(9)),]$Celltype_main_level <- "Smooth_muscle_cell"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(3)),]$Celltype_main_level <- "Fibroblast"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(17)),]$Celltype_main_level <- "Endothelial"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(12,10,8)),]$Celltype_main_level <- "Epithelial"
Seurat_object@meta.data[which(integrated_snn_res.0.4 %in% c(0,1,2,15,6,13)),]$Celltype_main_level <- "Eosinophils"
detach(Seurat_object@meta.data)
p <- DotPlot(Seurat_object, features = markers,group.by = "Celltype_main_level",)+
  RotatedAxis()+
  scale_x_discrete("")+
  scale_y_discrete("")
ggsave(filename = "./20231105_celltype/Other cluster annotation results.pdf",plot = p,width = 18,height = 8)

saveRDS(Seurat_object, "./20231105_celltype/Seurat_object.RDs")














