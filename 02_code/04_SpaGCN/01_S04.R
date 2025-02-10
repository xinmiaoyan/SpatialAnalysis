library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(cowplot)

setwd('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/03_seurat')

dir = '01_s04'
dir.create(dir)
setwd(dir)

data_dir <- '/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/02_spaceranger_count_s3/S-04-021004_without_align/S-04-021004/outs'
list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
obj = Load10X_Spatial(data.dir = data_dir)

#Data preprocessing
p1 <- VlnPlot(obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()+labs(x = NULL)+theme(axis.text.x = element_text(angle = 0, hjust = 1))
p2 <- SpatialFeaturePlot(obj, features = "nCount_Spatial") + theme(legend.position = "right")
p = plot_grid(p1, p2, ncol = 2, rel_widths = c(1/3, 2/3)) #wrap_plots(plotlist = list(p1, p2), layout_matrix = rbind(c(1), c(2,3,4)))
ggsave(plot = p, file = '01_qc.pdf', width = 7, height = 4)


#
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)

dir = '02_featurePlot'
dir.create(dir)
genes = c("PTPRC",'CD19','MS4A1','CD38','XBP1','JCHAIN',"CD3D",'CD4','FOXP3',
'CD8A','NCAM1','FCGR3A', 'CD68','CD163','MRC1','EPCAM','COL1A1','PECAM1','VWF')
##Gene expression visualization
for(i in genes){
p1 = SpatialFeaturePlot(obj, features = i)
# ggsave(plot = p, file = paste(dir,'/',i,'.pdf', sep=''), width = 6,height = 6)

# p1 <- SpatialFeaturePlot(obj, features = i, pt.size.factor = 1)
p2 <- SpatialFeaturePlot(obj, features = i, alpha = c(0.1, 1))
p = p1 + p2
ggsave(plot = p, file = paste(dir,'/',i,'.pdf', sep=''), width = 10,height = 5)

}


#Dimensionality reduction, clustering, and visualization
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 2)
p = p1 + p2
ggsave(plot = p, file = '02_clustering.pdf', width = 10,height = 5)

#Identification of Spatially Variable Features
de_markers <- FindMarkers(obj, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

cnts = obj@assays$Spatial@counts
write.table(cnts, "cnts.tsv", sep = "\t", row.names = T, quote = FALSE)

locs = obj@images$slice1@coordinates[c('imagerow','imagecol')]
colnames(locs) = c('y','x')
write.table(locs, "locs.tsv", sep = "\t", row.names = T, quote = FALSE)


