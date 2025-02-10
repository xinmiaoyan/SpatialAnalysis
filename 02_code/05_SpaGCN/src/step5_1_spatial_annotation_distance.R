#--------------------------------------------------------------
# filename : step5_spatial_cluster_score.R
# Date : 2023-05-26
# contributor : Yanshuo Chu
# function: step5_spatial_cluster_score
# R version: R/4.0.3
#--------------------------------------------------------------

print('<==== step5_spatial_cluster_score.R ====>')
rm(list=ls())

suppressMessages({
  library(optparse)
  library(tidyverse)
  library(SpatialClusterScore)
})

option_list = list(
  make_option(c("-c","--cellranger_out"),
                type = 'character',
                help = 'cellranger out',
                metavar = 'character'),
  make_option(c("-t","--annotation_table"),
              type = 'character',
              help = 'annotation table',
              metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

cor.matrix <- SpatialClusterScore::get_spatial_corr(
    cellranger_out = opt$cellranger_out,
    annotation_table = opt$annotation_table)

saveRDS(cor.matrix, file.path(dirname(opt$annotation_table), paste0('cor.matrix.rds')))
