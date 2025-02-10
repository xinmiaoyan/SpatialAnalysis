#--------------------------------------------------------------
# filename : step5_ha_plot.R
# Date : 2023-05-30
# contributor : Yanshuo Chu
# function: step5_ha_plot
# R version: R/4.0.3
#--------------------------------------------------------------

print("<==== step5_ha_plot.R ====>")
rm(list=ls())


suppressMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(RColorBrewer)
})

option_list = list(
    make_option(c("-t","--table"),
                type = 'character',
                help = 'spatial infor table',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

spatial_info_table <- read_tsv(opt$table)
## spatial_info_table <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/project26/result/PM3131/3_SpaGCN/HKU01/alpha_1_beta_49/p_0.4_nClusters_9/spatial_cluster_score_table.tsv")

color_table <- read_tsv("/rsrch3/scratch/genomic_med/ychu2/projects/project26/knowledge/private/pathology_annotation/pa_color.tsv")

spatial_info_table <- spatial_info_table %>%
  filter(!is.na(pathology_annotation))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unique(unlist(mapply(brewer.pal,
                                  qual_col_pals$maxcolors,
                                  rownames(qual_col_pals))))

gh <- ggplot(spatial_info_table) +
  geom_hex(aes(x = -y_array,
               y = -x_array,
               fill = pathology_annotation,
               color = pathology_annotation),
           lwd = 0.1,
           stat = "identity") +
  scale_fill_manual(values = color_table$Color,
                    breaks = color_table$Tissue) +
  scale_color_manual(values = color_table$Color,
                     breaks = color_table$Tissue) +
  theme_void()

nCluster = length(unique(spatial_info_table$cluster))
colorCluster = sample(col_vector, nCluster)
spatial_info_table$cluster <- as.character(spatial_info_table$cluster)
gc <- ggplot(spatial_info_table) +
  geom_hex(aes(x = -y_array,
               y = -x_array,
               fill = cluster,
               color = cluster),
           lwd = 0.1,
           stat = "identity") +
  scale_fill_manual(values = colorCluster) +
  scale_color_manual(values = colorCluster) +
  theme_void()

g <- cowplot::plot_grid(gh, gc, align = "hv", ncol = 2)

setwd(dirname(opt$table))
ggsave("histology_cluster_hex_plot.pdf", g, width = 12, height = 5)
