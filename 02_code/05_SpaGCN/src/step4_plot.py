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
import warnings

import scanpy as sc

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt


def plot_domain(adata_path, output_domain_plot_folder_path):
    adata = sc.read(adata_path)

    plot_color = [
        "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3", "#D1D1D1",
        "#6D1A9C", "#15821E", "#3A84E6", "#997273", "#787878", "#DB4C6C",
        "#9E7A7A", "#554236", "#AF5F3C", "#93796C", "#F9BD3F", "#DAB370",
        "#877F6C", "#268785"
    ]

    for domains in ["pred", "refined_pred"]:
        num_celltype = len(adata.obs[domains].unique())
        adata.uns[domains + "_colors"] = list(plot_color[:num_celltype])

        ax = sc.pl.scatter(adata,
                           alpha=1,
                           x="y_pixel",
                           y="x_pixel",
                           color=domains,
                           title=domains,
                           color_map=plot_color,
                           show=False,
                           size=100000 / adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()

        figure_path = "{0}/{1}.png".format(output_domain_plot_folder_path,
                                           domains)
        plt.savefig(figure_path, dpi=600)
        plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='cal adj')

    parser.add_argument('--adataPath', dest='adataPath', help='adata path')
    parser.add_argument('--outPath', dest='outPath', help='out path')

    args = parser.parse_args()

    print(args)

    plot_domain(adata_path=args.adataPath,
                output_domain_plot_folder_path=args.outPath)


if __name__ == '__main__':
    main()
