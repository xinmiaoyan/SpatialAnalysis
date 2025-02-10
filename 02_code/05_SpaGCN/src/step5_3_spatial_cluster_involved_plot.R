#--------------------------------------------------------------
# filename : step5_3_spatial_cluster_involved_plot.R
# Date : 2023-09-25
# contributor : Yanshuo Chu
# function: step5_3_spatial_cluster_involved_plot
# R version: R/4.0.3
#--------------------------------------------------------------

print("<==== step5_3_spatial_cluster_involved_plot.R ====>")
rm(list=ls())

suppressMessages({
  library(optparse)
  library(tidyverse)
  library(SpatialClusterScore)
})

option_list = list(
  make_option(c("-t","--annotation_table"),
              type = 'character',
              help = 'annotation table',
              metavar = 'character'),
  make_option(c("-c","--color_table"),
              type = 'character',
              help = 'color table',
              metavar = 'character'),
  make_option(c("-f","--pa_h_clusters_fill"),
              type = 'character',
              default = NULL,
              help = 'color table',
              metavar = 'character'),
  make_option(c("-b","--pa_h_clusters_border"),
              type = 'character',
              default = NULL,
              help = 'color table',
              metavar = 'character'),
  make_option(c("-m","--matrix"),
              type = 'character',
              help = 'filtered matrix h5 file path',
              metavar = 'character'),
  make_option(c("-n","--file_name"),
              type = 'character',
              default = NULL,
              help = 'color table',
              metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

st <- read_tsv(opt$annotation_table)

if(is.null(opt$pa_h_clusters_fill)) {
  SpatialClusterScore::plot_c_involved(spatial_info_table = opt$annotation_table,
                                       color_table = opt$color_table,
                                       output_file = opt$file_name,
                                       output_folder = dirname(opt$annotation_table))
}else{
  print(opt$pa_h_clusters_fill)
  print(opt$pa_h_clusters_border)

  if(is.null(opt$pa_h_clusters_fill)) {
    pa_h_c_f <- NULL
  }else{
    pa_h_c_f <- unlist(strsplit(opt$pa_h_clusters_fill, ","))
  }
  if(is.null(opt$pa_h_clusters_border)) {
    pa_h_c_b <- NULL
  }else{
    pa_h_c_b <- unlist(strsplit(opt$pa_h_clusters_border, ","))
  }

  print("pa_h_c_b")
  print(pa_h_c_b)
  print("pa_h_c_f")
  print(pa_h_c_f)
  print(all(pa_h_c_f %in% st$pathology_annotation))

  if(all(pa_h_c_f %in% st$pathology_annotation)) {
    print("file name  is")
    print(pa_h_c_f)
    SpatialClusterScore::plot_c_involved(spatial_info_table = opt$annotation_table,
                                         color_table = opt$color_table,
                                         exp_table = opt$matrix,
                                         output_folder = dirname(opt$annotation_table),
                                         output_file = opt$file_name,
                                         min_frac = 0.1,
                                         min_n = 10,
                                         highlight_PAs_fill = pa_h_c_f,
                                         highlight_PAs_border = pa_h_c_b)
  }
}
