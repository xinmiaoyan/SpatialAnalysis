#! /rsrch3/home/genomic_med/ychu2/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: step4_histology_annotation_ari.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
"""

import os

os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2, 40))
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score


def get_ari(adata_path, histology_annotation_path, sampleName):
    # adata_path = "/rsrch3/scratch/genomic_med/ychu2/projects/project26/result/PM3131/3_SpaGCN/HKU01/alpha_1_beta_49/p_0.5_nClusters_6/domain.h5ad"
    # histology_annotation_path = "/rsrch3/scratch/genomic_med/ychu2/projects/project26/knowledge/private/pathology_annotation/HKU01_Pathology_Annotations-type.csv"

    adata = sc.read(adata_path)
    histology_annotation_table = pd.read_csv(
        histology_annotation_path, header=0, na_filter=False, sep="\t"
    )

    adata.obs["Pathology_annotations"] = histology_annotation_table.loc[
        np.where(histology_annotation_table["Barcode"].isin(adata.obs.index)),
        "Pathology_annotations",
    ].values

    adata.obs["pathology_annotation_type"] = histology_annotation_table.loc[
        np.where(histology_annotation_table["Barcode"].isin(adata.obs.index)),
        "pathology_annotation_type",
    ].values

    spot_type_vector = adata.obs["Pathology_annotations"].to_list()
    cluster_vector = adata.obs["refined_pred"].to_list()

    alluvial_table = adata.obs[["Pathology_annotations", "pred", "refined_pred"]]
    alluvial_table.index.name = "spot"
    alluvial_table.to_csv(
        "{0}/{1}".format(os.path.dirname(adata_path), "alluvial.tsv"),
        sep="\t",
        index=True,
    )

    spatial_cluster_score_table = adata.obs[
        [
            "x_array",
            "y_array",
            "Pathology_annotations",
            "refined_pred",
            "pathology_annotation_type",
        ]
    ]
    spatial_cluster_score_table.index.name = "spot"
    spatial_cluster_score_table.rename(
        columns={
            "Pathology_annotations": "pathology_annotation",
            "refined_pred": "cluster",
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

    labelled_idx = [idx for idx, value in enumerate(spot_type_vector) if value != ""]

    ARI = adjusted_rand_score(
        [spot_type_vector[i] for i in labelled_idx],
        [cluster_vector[i] for i in labelled_idx],
    )

    outFilePath = "{0}/{1}".format(os.path.dirname(adata_path), "ARI.tsv")
    with open(outFilePath, "w") as outFile:
        outFile.write("{0}\t{2}\t{1}".format(outFilePath, ARI, sampleName))


def main():
    import argparse

    parser = argparse.ArgumentParser(description="")

    parser.add_argument("--adataPath", dest="adataPath", help="adata path")
    parser.add_argument(
        "--annotationTablePath",
        dest="annotationTablePath",
        help="pathology annotation table path",
    )
    parser.add_argument("--sampleName", dest="sampleName", help="sample name")

    args = parser.parse_args()

    print(args)

    get_ari(
        adata_path=args.adataPath,
        histology_annotation_path=args.annotationTablePath,
        sampleName=args.sampleName,
    )


if __name__ == "__main__":
    main()
