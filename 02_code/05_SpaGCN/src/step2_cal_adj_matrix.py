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

warnings.filterwarnings("ignore")


def cal_adj(adata_path, img_path, output_adj_file_path, alpha=1.0, beta=49.0):
    adata = sc.read(adata_path)
    img = cv2.imread(img_path)

    x_pixel = adata.obs["x_pixel"].tolist()
    y_pixel = adata.obs["y_pixel"].tolist()

    adj = spg.calculate_adj_matrix(x=x_pixel,
                                   y=y_pixel,
                                   x_pixel=x_pixel,
                                   y_pixel=y_pixel,
                                   image=img,
                                   beta=beta,
                                   alpha=alpha,
                                   histology=True)

    np.savetxt(output_adj_file_path, adj, delimiter=',')


def main():
    import argparse

    parser = argparse.ArgumentParser(description='cal adj')

    parser.add_argument('--adataPath', dest='adataPath', help='adata path')
    parser.add_argument('--imgPath', dest='imgPath', help='image path')
    parser.add_argument('--outPath', dest='outPath', help='out path')
    parser.add_argument('--alpha', dest='alpha', help='alpha')
    parser.add_argument('--beta', dest='beta', help='beta')

    args = parser.parse_args()

    print(args)

    cal_adj(adata_path=args.adataPath,
            img_path=args.imgPath,
            output_adj_file_path=args.outPath,
            alpha=float(args.alpha),
            beta=float(args.beta))


if __name__ == '__main__':
    main()
