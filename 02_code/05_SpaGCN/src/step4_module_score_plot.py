#! /rsrch3/home/genomic_med/ychu2/miniconda3/bin/python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: step4_module_score_plot.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
'''

import os

os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2, 40))
import scanpy as sc
import warnings

warnings.filterwarnings("ignore")
import pandas as pd
import matplotlib.pyplot as plt


def module_score_plot(adata_path, marker_table_path):
    adata = sc.read(adata_path)
    marker_table = pd.read_csv(marker_table_path,
                               header=0,
                               na_filter=False,
                               sep="\t")

    marker_table.columns.values[0] = "Marker"
    marker_table.columns.values[1] = "CellType"

    marker_set_name = os.path.splitext(os.path.basename(marker_table_path))[0]
    outDir = os.path.join(os.path.dirname(adata_path), marker_set_name)
    if not os.path.exists(outDir):
        os.mkdir(outDir, )

    for ct in marker_table['CellType'].unique():
        temp_marker = marker_table.loc[marker_table['CellType'] == ct,
                                       'Marker']
        sc.tl.score_genes(adata,
                          gene_list=temp_marker,
                          score_name="{}_score".format(ct))

        ax = sc.pl.scatter(adata,
                           alpha=1,
                           x="y_pixel",
                           y="x_pixel",
                           color="{}_score".format(ct),
                           title=ct,
                           show=False,
                           size=100000 / adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()

        figure_path = "{0}/{1}.pdf".format(outDir, ct)
        plt.savefig(figure_path)
        plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('--adataPath', dest='adataPath', help='adata path')
    parser.add_argument('--markerTablePath',
                        dest='markerTablePath',
                        help='marker table path')

    args = parser.parse_args()

    print(args)

    module_score_plot(adata_path=args.adataPath,
                      marker_table_path=args.markerTablePath)


if __name__ == '__main__':
    main()
