#--------------------------------------------------------------
# filename : step5_post_corr.R
# Date : 2023-05-30
# contributor : Yanshuo Chu
# function: step5_post_corr
# R version: R/4.2.1
#--------------------------------------------------------------

print("<==== step5_post_corr.R ====>")
rm(list=ls())

suppressMessages({
  library(optparse)
  library(Seurat)
  library(tidyverse)
  library(pheatmap)
})

option_list = list(
    make_option(c("-d","--data"),
                type = 'character',
                help = 'out folder',
                metavar = 'character'),
    make_option(c("-a","--annotation"),
                type = 'character',
                help = 'annotation',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

seurat_obj <- Seurat::Load10X_Spatial(opt$data)
Annotation_T <- read_tsv(opt$annotation)

seurat_obj <- subset(seurat_obj, cells = Annotation_T$spot)
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial")
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, selection.method = "moransi")
svf <- Seurat::SpatiallyVariableFeatures(seurat_obj, selection.method = "moransi")
seurat_obj_svf <- subset(seurat_obj, features = svf)

seurat_obj_svf$pathology_annotation <-
  Annotation_T$pathology_annotation[
                 match(Cells(seurat_obj_svf), Annotation_T$spot)]
seurat_obj_svf$cluster <-
  Annotation_T$cluster[match(Cells(seurat_obj_svf), Annotation_T$spot)]

Idents(seurat_obj_svf) <- seurat_obj_svf$pathology_annotation

avg.exp <- AverageExpression(seurat_obj_svf, slot = "data")$Spatial
avg.exp <- avg.exp[!is.na(colSums(avg.exp)),]
cor.exp <- as.matrix(cor(avg.exp, method = "spearman"))
pheatmap::pheatmap(cor.exp,
                   cluster_rows = T,
                   cluster_cols = T,
                   filename = file.path(dirname(opt$annotation),
                                        "pathology_annotation_spatial_corr.pdf"))

saveRDS(cor.exp, file.path(dirname(opt$annotation), paste0('ha_cor_exp.rds')))

## min_v <- min(cor.exp)
## max_v <- max(cor.exp)
## cor.exp <- (max_v - cor.exp) / (max_v - min_v)
## saveRDS(cor.exp, file.path(getwd(), paste0('cor.exp', "_", Sys.Date(), '.rds')))
## c <- 0
## ct <- table(seurat_obj@meta.data$pathology_annotation, seurat_obj@meta.data$cluster)
## for(i in 1:dim(ct)[1]) {
##   for(j in 1:dim(ct)[1]) {
##     c <- c + sum(ct[i,] * ct[j,] * cor.exp[i,j]) / 2
##   }
## }

Idents(seurat_obj_svf) <- seurat_obj_svf$cluster
avg.exp <- AverageExpression(seurat_obj_svf, slot = "data")$Spatial
avg.exp <- avg.exp[!is.na(colSums(avg.exp)),]
cor.exp <- as.matrix(cor(avg.exp, method = "spearman"))
pheatmap::pheatmap(cor.exp,
                   cluster_rows = T,
                   cluster_cols = T,
                   filename = file.path(dirname(opt$annotation),
                                        "cluster_spatial_corr.pdf"))

saveRDS(cor.exp, file.path(dirname(opt$annotation), paste0('cluster_cor_exp.rds')))

