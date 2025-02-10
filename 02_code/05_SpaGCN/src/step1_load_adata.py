# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: load.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
'''

import warnings

import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

print("pandas version {}".format(pd.__version__))


def load_10x_to_adata(filtered_feature_bc_matrix_path, tissue_position_path,
                      output_adata_file_path):
    # adata = sc.read_10x_h5(
    #     "/rsrch3/scratch/genomic_med/ychu2/projects/project26/result/PM3131/2_spaceranger_count/HKU01/outs/filtered_feature_bc_matrix.h5"
    # )
    # spatial = pd.read_csv(
    #     "/rsrch3/scratch/genomic_med/ychu2/projects/project26/result/PM3131/2_spaceranger_count/HKU01/outs/spatial/tissue_positions.csv",
    # header=0,
    # na_filter=False,
    # index_col=0)
    adata = sc.read_10x_h5(filtered_feature_bc_matrix_path)
    spatial = pd.read_csv(tissue_position_path,
                          header=0,
                          na_filter=False,
                          index_col=0)

    for i in range(5):
        spatial[spatial.columns[i]] = \
            pd.to_numeric(spatial[spatial.columns[i]])

    adata.obs["x1"] = spatial[spatial.columns[0]]
    adata.obs["x2"] = spatial[spatial.columns[1]]
    adata.obs["x3"] = spatial[spatial.columns[2]]
    adata.obs["x4"] = spatial[spatial.columns[3]]
    adata.obs["x5"] = spatial[spatial.columns[4]]
    adata.obs["x_array"] = adata.obs["x2"]
    adata.obs["y_array"] = adata.obs["x3"]
    adata.obs["x_pixel"] = adata.obs["x4"]
    adata.obs["y_pixel"] = adata.obs["x5"]
    adata = adata[adata.obs["x1"] == 1]
    adata.var_names = [i.upper() for i in list(adata.var_names)]
    adata.var["genename"] = adata.var.index.astype("str")
    adata.write_h5ad(output_adata_file_path)


def main():
    import argparse

    parser = argparse.ArgumentParser(description='load to adata')

    parser.add_argument('--matrixPath', dest='matrixPath', help='matrix path')
    parser.add_argument('--tissuePositionPath',
                        dest='tissuePositionPath',
                        help='tissue position path')
    parser.add_argument('--outPath', dest='outPath', help='out path')

    args = parser.parse_args()
    print(args)

    load_10x_to_adata(filtered_feature_bc_matrix_path=args.matrixPath,
                      tissue_position_path=args.tissuePositionPath,
                      output_adata_file_path=args.outPath)


if __name__ == '__main__':
    main()
