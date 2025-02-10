# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: step4_DEG.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
"""

import argparse
import os

import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def get_DEG_by_group(adata, out_path, groupby="clusters"):
    """get DEG by group"""
    sc.tl.rank_genes_groups(adata=adata, groupby=groupby, method="t-test")
    DEGs_group = adata.uns["rank_genes_groups"]
    groups = DEGs_group["names"].dtype.names

    df = None
    for group in groups:
        temp_df = pd.DataFrame(
            {
                key: DEGs_group[key][group]
                for key in ["names", "logfoldchanges", "pvals", "pvals_adj"]
            }
        )
        temp_df["cluster"] = group
        if df is None:
            df = temp_df
        else:
            df = pd.concat([df, temp_df])

    gene_pattern = [
        "MALAT1",
        "^MT-",
        "^RPL",
        "^RPS",
        "^LOC(0-9)",
        "^TR(A|B|G|D)V",
        "^MTRNR",
    ]
    dft = df.names.str.contains(
        "(" + "|".join(gene_pattern) + ")", regex=True, na=False
    )
    df = df.loc[np.invert(dft), :]
    df.to_csv(os.path.join(out_path, "DEGs_groupby_{}.tsv".format(groupby)), sep="\t")

    sc.pl.rank_genes_groups_heatmap(
        adata=adata,
        n_genes=20,
        groupby=groupby,
        show_gene_labels=True,
        swap_axes=True,
        cmap="gnuplot2",
    )
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, "DEGs_heatmap_groupby_{}.pdf".format(groupby)))

    sc.pl.rank_genes_groups_tracksplot(adata, n_genes=20, groupby=groupby)
    plt.tight_layout()
    plt.savefig(
        os.path.join(out_path, "DEGs_tracksplot_groupby_{}.pdf".format(groupby))
    )

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=20, groupby=groupby, cmap="gnuplot2")
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, "DEGs_dotplot_groupby_{}.pdf".format(groupby)))

    sc.pl.rank_genes_groups_stacked_violin(
        adata, n_genes=20, groupby=groupby, swap_axes=True
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(out_path, "DEGs_stackedviolin_groupby_{}.pdf".format(groupby))
    )

    sc.pl.correlation_matrix(adata=adata, groupby=groupby)
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, "correlation_groupby_{}.pdf".format(groupby)))


def run_rank_gene(adata_path, annotation_table_path, out_path):
    """generate DEGs"""
    adata = sc.read(adata_path)

    annotation_table = pd.read_csv(annotation_table_path, header=0, na_filter=False)
    annotation_table["Pathology_annotations"] = (
        annotation_table["Pathology_annotations"].str.strip().replace("", "Unlabelled")
    )
    tempObs = adata.obs.merge(
        annotation_table, left_index=True, right_on="Barcode", how="left"
    )
    tempObs.set_index("Barcode", inplace=True)
    adata.obs = tempObs

    adata = run_umap(adata, annotation_table_path, out_path)

    for groupby in ["refined_pred", "Pathology_annotations", "clusters"]:
        adata.obs[groupby] = adata.obs[groupby].astype(str)
        temp_out_path = os.path.join(out_path, "DEGs_groupby_{}".format(groupby))
        if not os.path.exists(temp_out_path):
            os.mkdir(temp_out_path)
        get_DEG_by_group(adata, temp_out_path, groupby=groupby)


def run_umap(adata, annotation_table_path, out_path):
    """run umap"""
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="clusters")


    annotation_table = pd.read_csv(annotation_table_path,
                                   header=0,
                                   na_filter=False)
    annotation_table['Pathology_annotations'] = annotation_table[
        'Pathology_annotations'].str.strip().replace('', 'Unlabelled')
    unique_PAs = annotation_table["Pathology_annotations"].unique()

    color_df = pd.read_csv("/rsrch3/scratch/genomic_med/ychu2/projects/project26/knowledge/private/pathology_annotation/pa_color.tsv", sep='\t')
    color_values = [
        color_df.loc[color_df['Tissue'] == t, 'Color'].iloc[0]
        for t in unique_PAs
    ]
    color_dict1 = dict(zip(unique_PAs, color_values))
    color_dict2 = dict(zip([str(i) for i in range(20)], sc.pl.palettes.default_20))
    color_dict3 = dict(zip([i for i in range(20)], sc.pl.palettes.default_20))
    color_dict1.update(color_dict2)
    color_dict1.update(color_dict3)
    print(color_dict1)

    plt.rcParams["figure.figsize"] = (4, 4)
    sc.pl.umap(
        adata,
        color=["refined_pred", "Pathology_annotations", "clusters"],
        palette=color_dict1,
        wspace=0.8,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, "UMAP.pdf"))

    return adata


def main():
    parser = argparse.ArgumentParser(description="clustering")
    parser.add_argument("--adataPath", dest="adata_path", help="adata object")
    parser.add_argument("--outPath", dest="out_path", help="out path")
    parser.add_argument(
        "--annotationTablePath",
        dest="annotation_table_path",
        help="annotation table path",
    )

    args = parser.parse_args()
    print(args)

    run_rank_gene(
        adata_path=args.adata_path,
        annotation_table_path=args.annotation_table_path,
        out_path=args.out_path,
    )


if __name__ == "__main__":
    main()
