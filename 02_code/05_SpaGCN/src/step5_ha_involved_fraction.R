#--------------------------------------------------------------
# filename : step5_ha_involved_fraction.R
# Date : 2023-09-25
# contributor : Yanshuo Chu
# function: step5_ha_involved_fraction
# R version: R/4.0.3
#--------------------------------------------------------------

print("<==== step5_ha_involved_fraction.R ====>")
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
  make_option(c("-n","--min_n"),
              type = 'integer',
              default = 10,
              help = 'minimum number of c involved',
              metavar = 'integer'),
  make_option(c("-f","--min_frac"),
              type = 'double',
              default = 0.1,
              help = 'minimum fraction of c involved',
              metavar = 'double')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

## at <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/project27/result/kevin_data/GARP/SpaGCN/IX_56689/alpha_1_beta_49/p_0.5_nClusters_17/spatial_cluster_score_table.tsv")
at <- read_tsv(opt$annotation_table)

pas <- unique(at$pathology_annotation)

frac_t <- c()

for(pa in pas) {
  c_involved <- at %>%
    filter(pathology_annotation == pa) %>%
    group_by(cluster) %>%
    count() %>%
    ungroup() %>% 
    mutate(frac = n / sum(n)) %>%
    filter(n > opt$min_n | frac > opt$min_frac) %>% 
    pull(cluster) %>%
    unique()

  c_involved_frac <- at %>%
    filter(cluster %in% c_involved |
           pathology_annotation == pa) %>%
    nrow() / at %>%
    filter(pathology_annotation == pa) %>%
    nrow()
  
  temp_t <- tibble(path = opt$annotation_table,
                   min_n = opt$min_n,
                   min_frac = opt$min_frac,
                   pa = pa,
                   c_involved_frac = c_involved_frac)

  frac_t <- bind_rows(frac_t, temp_t)
}


combined_list <- list(mic = c("Micrometastasis Stroma", "Micrometastasis"),
                      tumor = c("Stroma rich tumor bed", "Tumor", "Immune rich stroma tumor bed"))
for (i in seq_along(combined_list)) {
  temp_name <- names(combined_list[i])
  temp_list <- combined_list[[i]]
  o_l <- intersect(temp_list, pas)
  if (length(o_l) > 1) {
    c_involved <- at %>%
      filter(pathology_annotation %in% o_l) %>%
      group_by(cluster) %>%
      count() %>%
      ungroup() %>% 
      mutate(frac = n / sum(n)) %>%
      filter(n > opt$min_n | frac > opt$min_frac) %>% 
      pull(cluster) %>%
      unique()

    c_involved_frac <- at %>%
      filter(cluster %in% c_involved |
             pathology_annotation %in% o_l) %>%
      nrow() / at %>%
      filter(pathology_annotation %in% o_l) %>%
      nrow()
    
    temp_t <- tibble(path = opt$annotation_table,
                     min_n = opt$min_n,
                     min_frac = opt$min_frac,
                     pa = paste0(o_l, collapse = "-"),
                     c_involved_frac = c_involved_frac)

    frac_t <- bind_rows(frac_t, temp_t)
  }
}

write_tsv(frac_t, file.path(dirname(opt$annotation_table), paste0('c_involved_frac.tsv')))

