# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: step2_cal_adj_matrix.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
'''
import os

os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2, 40))

import random
import warnings

import numpy as np
import scanpy as sc
import SpaGCN as spg
import torch

warnings.filterwarnings("ignore")


def detect_domain(adata_path,
                  adj_path,
                  output_adj_file_path,
                  p=0.5,
                  n_clusters=7,
                  seq_type="Visium"):
    adata = sc.read(adata_path)
    adj = np.loadtxt(adj_path, delimiter=',')

    adata.var_names_make_unique()
    spg.prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)

    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    r_seed = t_seed = n_seed = 100

    res = spg.search_res(adata,
                         adj,
                         l,
                         n_clusters,
                         start=0.7,
                         step=0.1,
                         tol=5e-3,
                         lr=0.05,
                         max_epochs=20,
                         r_seed=r_seed,
                         t_seed=t_seed,
                         n_seed=n_seed)

    clf = spg.SpaGCN()
    clf.set_l(l)

    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)

    #Run
    clf.train(adata,
              adj,
              init_spa=True,
              init="louvain",
              res=res,
              tol=5e-3,
              lr=0.05,
              max_epochs=200)

    y_pred, prob = clf.predict()

    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype('category')

    #Do cluster refinement(optional)
    #shape="hexagon" for Visium data, "square" for ST data.

    shape_str = "hexagon"
    if seq_type != "Visium":
        shape_str = "square"

    x_array = adata.obs["x_array"].tolist()
    y_array = adata.obs["y_array"].tolist()
    adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)

    refined_pred = spg.refine(sample_id=adata.obs.index.tolist(),
                              pred=adata.obs["pred"].tolist(),
                              dis=adj_2d,
                              shape=shape_str)

    adata.obs["refined_pred"] = refined_pred
    adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')

    #Save results
    adata.write_h5ad(output_adj_file_path)


def main():
    import argparse

    parser = argparse.ArgumentParser(description='cal adj')

    parser.add_argument('--adataPath', dest='adataPath', help='adata path')
    parser.add_argument('--adjPath', dest='adjPath', help='adj path')
    parser.add_argument('--outPath', dest='outPath', help='out path')
    parser.add_argument(
        '--p',
        dest='p',
        help='percentage of total expression contribution by neighborhoods, 0-1'
    )
    parser.add_argument('--nClusters',
                        dest='nClusters',
                        help='number of cluster')

    args = parser.parse_args()

    print(args)

    detect_domain(adata_path=args.adataPath,
                  adj_path=args.adjPath,
                  output_adj_file_path=args.outPath,
                  p=float(args.p),
                  n_clusters=int(args.nClusters))


if __name__ == '__main__':
    main()
