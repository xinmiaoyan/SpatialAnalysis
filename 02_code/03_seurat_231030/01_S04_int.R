

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(cowplot)

setwd('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/03_seurat')

dir = '01_S04'
dir.create(dir)
setwd(dir)

data_dir <- '/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/02_sapceranger_count_S04_S11_aglin/S04_align_Spaceranger'
list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
obj = Load10X_Spatial(data.dir = data_dir)

# add meta data
df1 = read.csv('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/03_seurat/01_S04/00_data/S04_ANNOTATION.csv')
df2 = read.csv('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/03_seurat/01_S04/00_data/S04_HEV.csv')
df3 = read.csv('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/16_SP_TLS_jianfeng/03_result/03_seurat/01_S04/00_data/S04_Neutrophils.csv')

table(df1$ANNOTATION)
table(df2$HEV)
table(df3$Neutrophils)
# > table(df1$ANNOTATION)

#                Diffuse TILs   Neutrophils Normal stroma        Stroma
#          1579            73            35          3250          1999
#           TLS         Tumor    Urothelium       Vessels
#            88          1005           168           832
# > table(df2$HEV)

#       HEV
# 9018   11
# > table(df3$Neutrophils)

#             Neutrophils
#        9004          25

# d1 = as.list(df1 %>% filter(ANNOTATION=='Neutrophils'))$Barcode
# d2 = as.list(df3 %>% filter(Neutrophils=='Neutrophils'))$Barcode
# con = intersect(d1,d2)
# con = setdiff(d1,d2)

identical(rownames(obj@meta.data), colnames(obj))
identical(df1$Barcode, colnames(obj)) #True
obj$Annotation = df1$ANNOTATION

#false
# df1$Barcode = factor(df1$Barcode, levels = colnames(obj))



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

for(dim in seq(20,50,10)){
    for(res in seq(0.3,1,0.2)){
        dir = paste('dim_',dim,'res_'res)
        dir.create(dir)
        setwd(dir)

        obj <- FindNeighbors(obj, reduction = "pca", dims = 1:dim)
        obj <- FindClusters(obj, resolution = res,verbose = FALSE)
        obj <- RunUMAP(obj, reduction = "pca",  dims = 1:dim)        

        p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
        p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 2)
        p3 <- SpatialDimPlot(obj, group.by = 'Annotation')
        p = p1 + p2 + p3
        ggsave(plot = p, file = paste('02_clustering_',dim,'_',res,'.pdf',sep=''), width = 15,height = 5)

        p = SpatialDimPlot(obj, cells.highlight = CellsByIdentities(obj, idents = levels(unique(Idents(obj)))), facet.highlight = TRUE, ncol = 3)
        ggsave(plot = p, file = paste('03_clustering_',dim,'_',res,'.pdf',sep=''), width = 15,height = 5)

        obj <- FindSpatiallyVariableFeatures(obj, assay = "SCT", features = VariableFeatures(obj)[1:1000],selection.method = "moransi")
        
        sp_markers = SpatiallyVariableFeatures(obj, selection.method = "moransi")
        write_csv(sp_markers, '04_sp_markers.csv')
        top.features <- head(SpatiallyVariableFeatures(obj, selection.method = "moransi"), 6)
        p = SpatialFeaturePlot(obj, features = top.features, ncol = 3, alpha = c(0.1, 1))
        ggsave(plot = p, file = paste('04_sp_marker_',dim,'_',res,'.pdf',sep=''), width = 15,height = 5)

        markers <- FindAllMarkers(obj)
        markers = markers%>%filter(p_val_adj<0.05)
        markers = arrange(markers, cluster, desc(avg_log2FC), desc(p_val_adj))
        write_csv(markers, '05_sc_markers.csv')

        setwd('../')

    }
}




# #Identification of Spatially Variable Features
# de_markers <- FindMarkers(obj, ident.1 = 5, ident.2 = 6)
# SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)


# cnts = obj@assays$Spatial@counts
# write.table(cnts, "cnts.tsv", sep = "\t", row.names = T, quote = FALSE)

# locs = obj@images$slice1@coordinates[c('imagerow','imagecol')]
# colnames(locs) = c('y','x')
# write.table(locs, "locs.tsv", sep = "\t", row.names = T, quote = FALSE)

















