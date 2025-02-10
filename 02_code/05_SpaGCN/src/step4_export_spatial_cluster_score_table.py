# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: step4_histology_annotation_ari.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
"""

import argparse
import os

os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2, 40))
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score


def export_table(adata_path, histology_annotation_path):
    # adata_path = "/rsrch3/scratch/genomic_med/ychu2/projects/project27/result/kevin_data/GARP/SpaGCN/IX_56688/alpha_1.5_beta_70/p_0.6_nClusters_9/domain.h5ad"
    # histology_annotation_path = "/rsrch3/scratch/genomic_med/ychu2/projects/project27/data/kevin_data/GARP/pa/pa_IX_56688.csv"

    adata = sc.read(adata_path)
    histology_annotation_table = pd.read_csv(
        histology_annotation_path, header=0, na_filter=False, sep=","
    )
    histology_annotation_table.columns = ["Barcode", "Pathology_annotations"]

    if 'Pathology_annotations' in adata.obs:
        adata.obs.drop('Pathology_annotations', axis = 1)

    adata.obs = pd.merge(
        adata.obs,
        histology_annotation_table,
        left_index=True,
        right_on="Barcode",
        how="left",
    )
    adata.obs.set_index("Barcode", inplace=True)

    spatial_cluster_score_table = adata.obs[
        [
            "x_array",
            "y_array",
            "Pathology_annotations",
            "refined_pred",
        ]
    ]
    spatial_cluster_score_table.index.name = "spot"
    spatial_cluster_score_table.rename(
        columns={
            "Pathology_annotations": "pathology_annotation",
            "refined_pred": "cluster",
            "spot": "Barcode",
        },
        inplace=True,
    )

    spatial_cluster_score_table.to_csv(
        "{0}/{1}".format(
            os.path.dirname(adata_path), "spatial_cluster_score_table.tsv"
        ),
        sep="\t",
        index=True,
    )


def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("--adataPath", dest="adataPath", help="adata path")
    parser.add_argument(
        "--annotationTablePath",
        dest="annotationTablePath",
        help="pathology annotation table path",
    )

    args = parser.parse_args()

    print(args)

    export_table(
        adata_path=args.adataPath,
        histology_annotation_path=args.annotationTablePath,
    )


if __name__ == "__main__":
    main()
