# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: step2_cal_adj_matrix.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
"""

import argparse
import os

os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2, 40))
import csv
import math
import random
import re
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import SpaGCN as spg
import torch
from scipy.sparse import issparse

warnings.filterwarnings("ignore")
import cv2
import matplotlib.colors as clr
import matplotlib.pyplot as plt


def detect_SVG(raw_adata_path, adata_path, output_folder_path, target=0):
    raw = sc.read(raw_adata_path)
    adata = sc.read(adata_path)

    raw.var_names_make_unique()
    raw.obs["pred"] = adata.obs["pred"].astype("category")
    raw.obs["x_array"] = raw.obs["x2"]
    raw.obs["y_array"] = raw.obs["x3"]
    raw.obs["x_pixel"] = raw.obs["x4"]
    raw.obs["y_pixel"] = raw.obs["x5"]
    # Convert sparse matrix to non-sparse
    raw.X = raw.X.A if issparse(raw.X) else raw.X
    # raw.raw=raw
    sc.pp.log1p(raw)

    # Set filtering criterials
    min_in_group_fraction = 0.8
    min_in_out_group_ratio = 1
    min_fold_change = 1.5

    x_array = adata.obs["x_array"].tolist()
    y_array = adata.obs["y_array"].tolist()

    # Search radius such that each spot in the target domain has approximately 10 neighbors on average
    adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    start, end = np.quantile(adj_2d[adj_2d != 0], q=0.001), np.quantile(
        adj_2d[adj_2d != 0], q=0.1
    )
    r = spg.search_radius(
        target_cluster=target,
        cell_id=adata.obs.index.tolist(),
        x=x_array,
        y=y_array,
        pred=adata.obs["pred"].tolist(),
        start=start,
        end=end,
        num_min=10,
        num_max=14,
        max_run=100,
    )
    # Detect neighboring domains
    nbr_domians = spg.find_neighbor_clusters(
        target_cluster=target,
        cell_id=raw.obs.index.tolist(),
        x=raw.obs["x_array"].tolist(),
        y=raw.obs["y_array"].tolist(),
        pred=raw.obs["pred"].tolist(),
        radius=r,
        ratio=1 / 2,
    )

    nbr_domians = nbr_domians[0:3]
    de_genes_info = spg.rank_genes_groups(
        input_adata=raw,
        target_cluster=target,
        nbr_list=nbr_domians,
        label_col="pred",
        adj_nbr=True,
        log=True,
    )
    # Filter genes
    de_genes_info = de_genes_info[(de_genes_info["pvals_adj"] < 0.05)]
    filtered_info = de_genes_info
    filtered_info = filtered_info[
        (filtered_info["pvals_adj"] < 0.05)
        & (filtered_info["in_out_group_ratio"] > min_in_out_group_ratio)
        & (filtered_info["in_group_fraction"] > min_in_group_fraction)
        & (filtered_info["fold_change"] > min_fold_change)
    ]
    filtered_info = filtered_info.sort_values(by="in_group_fraction", ascending=False)
    filtered_info["target_dmain"] = target
    filtered_info["neighbors"] = str(nbr_domians)
    print("SVGs for domain ", str(target), ":", filtered_info["genes"].tolist())

    file_path = "{0}/domain_{1}_SVG.tsv".format(output_folder_path, target)
    filtered_info.to_csv(file_path, sep="\t")

    # Plot refinedspatial domains
    color_self = clr.LinearSegmentedColormap.from_list(
        "pink_green", ["#3AB370", "#EAE7CC", "#FD1593"], N=256
    )
    # for g in filtered_info["genes"].tolist():
    #     raw.obs["exp"] = raw.X[:, raw.var.index == g]
    #     ax = sc.pl.scatter(raw,
    #                        alpha=1,
    #                        x="y_pixel",
    #                        y="x_pixel",
    #                        color="exp",
    #                        title=g,
    #                        color_map=color_self,
    #                        show=False,
    #                        size=100000 / raw.shape[0])
    #     ax.set_aspect('equal', 'box')
    #     ax.axes.invert_yaxis()
    #     plt.savefig(output_folder_path + "/" + g + ".png", dpi=600)
    #     plt.close()


def main():
    parser = argparse.ArgumentParser(description="cal adj")

    parser.add_argument("--rawAdataPath", dest="rawAdataPath", help="raw adata path")
    parser.add_argument("--adataPath", dest="adataPath", help="adata path")
    parser.add_argument("--outDirPath", dest="outDirPath", help="out dir path")
    parser.add_argument("--target", dest="target", help="target domain")

    args = parser.parse_args()

    print(args)

    detect_SVG(
        raw_adata_path=args.rawAdataPath,
        adata_path=args.adataPath,
        output_folder_path=args.outDirPath,
        target=int(args.target),
    )


if __name__ == "__main__":
    main()
