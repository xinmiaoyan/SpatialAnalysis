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
  library(Seurat)
  library(SpatialClusterScore)
})


option_list = list(
  make_option(c("-d","--distance"),
              type = 'character',
              help = 'distance rds',
              metavar = 'character'),
  make_option(c("-a","--alpha"),
              type = 'double',
              default = 1.0,
              help = 'alpha',
              metavar = 'double'),
  make_option(c("-b","--beta"),
              type = 'double',
              default = 0.0,
              help = 'beta',
              metavar = 'double'),
  make_option(c("-m","--matrix"),
              type = 'character',
              help = 'filtered matrix h5 file path',
              metavar = 'character'),
  make_option(c("-t","--annotation_table"),
              type = 'character',
              help = 'annotation table',
              metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

cor.matrix <- readRDS(opt$distance)

score <- SpatialClusterScore::SpatialClusterScore(
                                D = cor.matrix,
                                spatial_info_table = opt$annotation_table,
                                alpha = opt$alpha,
                                beta = opt$beta)

cluster_rogue <- SpatialClusterScore::SpatialRogueScore(
                                exp_table = opt$matrix,
                                spatial_info_table = opt$annotation_table)

ha_rogue <- SpatialClusterScore::SpatialRogueScore_ha(
                                exp_table = opt$matrix,
                                spatial_info_table = opt$annotation_table)


tempT <- tibble(
  path = dirname(opt$annotation_table),
  SCS = score,
  cluster_ROGUE = cluster_rogue,
  ha_ROGUE = ha_rogue,
  alpha = opt$alpha,
  beta = opt$beta)

write_tsv(tempT, file.path(dirname(opt$distance),
                           paste0('SCS_alpha_', opt$alpha, '_beta_', opt$beta, '_.tsv')))
