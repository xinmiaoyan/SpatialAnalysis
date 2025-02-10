#! /rsrch3/home/genomic_med/ychu2/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: step4_plot_gene.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
"""

import os

os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2, 40))
import warnings

import scanpy as sc

warnings.filterwarnings("ignore")
import math

import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial.distance import cdist


def feature_plot(adata_path, marker_table_path, out_folder_path):
    adata = sc.read(adata_path)

    if not os.path.exists(out_folder_path):
        os.mkdir(out_folder_path)

    marker_table = pd.read_csv(marker_table_path, header=0, na_filter=False, sep="\t")
    marker_table.columns.values[0] = "Marker"

    for marker in marker_table["Marker"].unique():
        ax = sc.pl.scatter(
            adata,
            alpha=1,
            x="y_pixel",
            y="x_pixel",
            color=marker,
            title=ct,
            show=False,
            size=100000 / adata.shape[0],
        )
        ax.set_aspect("equal", "box")
        ax.axes.invert_yaxis()

        figure_path = "{0}/{1}.pdf".format(out_folder_path, marker)
        plt.savefig(figure_path)
        plt.close()

    pass


def feature_plot_hex(adata_path, marker_table_path, out_folder_path):
    # out_folder_path = "."
    # marker_table_path = "/rsrch3/scratch/genomic_med/ychu2/projects/project26/knowledge/private/GastricGeneOfInterest.tsv"
    # adata_path = "/rsrch3/scratch/genomic_med/ychu2/projects/project26/result/PM3131/3_SpaGCN/HKU01/alpha_1_beta_49/p_0.6_nClusters_10/domain.h5ad"
    adata = sc.read(adata_path)
    if not os.path.exists(out_folder_path):
        os.mkdir(out_folder_path)
    marker_table = pd.read_csv(marker_table_path, header=0, na_filter=False, sep="\t")
    marker_table.columns.values[0] = "Marker"

    for marker in marker_table["Marker"].unique():
        if marker in adata.var_names:
            plt.figure(figsize=(8, 7), linewidth=0.0)
            hb = plt.hexbin(
                max(adata.obs["y_array"]) - adata.obs["y_array"],
                max(adata.obs["x_array"]) - adata.obs["x_array"],
                C=adata[adata.obs.index, marker].X.toarray().flatten(),
                gridsize=(
                    math.floor(
                        (max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2
                    ),
                    math.floor(
                        (max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2
                    ),
                ),
                edgecolors="face",
                linewidth=0.1,
                cmap="coolwarm",
            )
            plt.axis("off")
            plt.gca().set_frame_on(False)
            cb = plt.colorbar(hb, shrink=0.6)
            cb.set_label("log1pEXP")
            figure_path = "{0}/{1}.pdf".format(out_folder_path, marker)
            # figure_path = "test.pdf"
            plt.title(marker)
            plt.tight_layout()
            plt.savefig(figure_path)
            plt.close()

            pass
        pass
    pass


def feature_plot_hex_ha(
    adata_path,
    marker_table_path,
    annotation_table_path,
    color_table_path,
    out_folder_path,
    marker_index=0,
):
    color_df = pd.read_csv(color_table_path, sep="\t")

    if isinstance(adata_path, str):
        adata = sc.read(adata_path)
    else:
        adata = adata_path

    if not os.path.exists(out_folder_path):
        os.mkdir(out_folder_path)
    marker_table = pd.read_csv(marker_table_path, header=0, na_filter=False, sep="\t")
    marker_table.columns.values[marker_index] = "Marker"

    annotation_table = pd.read_csv(annotation_table_path, header=0, na_filter=False)

    spot_array = np.array(
        [
            [x, y]
            for x, y in zip(
                (max(adata.obs["y_array"]) - adata.obs["y_array"]),
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            )
        ]
    )

    for marker in marker_table["Marker"].unique():
        if marker in adata.var_names:
            plt.figure(figsize=(8, 7), linewidth=0.0)
            hb = plt.hexbin(
                max(adata.obs["y_array"]) - adata.obs["y_array"],
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
                C=adata[adata.obs.index, marker].X.toarray().flatten(),
                gridsize=(
                    math.floor(
                        (max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2
                    ),
                    math.floor(
                        (max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2
                    ),
                ),
                edgecolors="face",
                linewidth=0.1,
                cmap="gnuplot2",
            )
            plt.axis("off")
            plt.gca().set_frame_on(False)
            cb = plt.colorbar(hb, shrink=0.6)
            cb.set_label("log1pEXP")
            plt.title(marker)
            plt.tight_layout()
            verts = hb.get_offsets()
            unique_y = np.unique(np.array([item[1] for item in verts]))
            diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
            unique_x = np.unique(np.array([item[0] for item in verts]))
            diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
            index = get_nearest_index(verts, spot_array)
            pa_indexed = adata.obs.index[index]
            index_df = pd.Index(pa_indexed)
            pa = (
                annotation_table.loc[annotation_table["Barcode"].isin(pa_indexed)]
                .set_index("Barcode")
                .reindex(index_df)["Pathology_annotations"]
                .str.strip()
            )
            for tpa in np.unique(pa):
                if tpa != "":
                    tcolor = color_df.loc[color_df["Tissue"] == tpa, "Color"].values[0]
                    plot_outline(
                        verts=verts,
                        pa=pa,
                        tpa=tpa,
                        plt=plt,
                        diff_y=diff_y,
                        diff_x=diff_x,
                        lw=0.5,
                        beta=0.88,
                        color=tcolor,
                    )
            figure_path = "{0}/{1}.pdf".format(out_folder_path, marker)
            plt.savefig(figure_path)
            plt.close()


def plot_outline(verts, pa, tpa, plt, diff_y, diff_x, lw, color="white", beta=0.88):
    nv = verts[pa == tpa]
    lines = np.concatenate(
        [hexLines(a=diff_y * 2 / 3, i=diff_x, off=off) for off in nv]
    )
    directions = np.concatenate([getDirections(off=off) for off in nv])
    uls, c = np.unique(lines.round(2), axis=0, return_counts=True)
    reshaped_array1 = uls[c == 1][:, np.newaxis]
    uls1_index = np.where(np.all(lines.round(2) == reshaped_array1, axis=(2, 3)))[1]
    t_d = directions[uls1_index]
    t_d = t_d[np.lexsort((t_d[:, 1], t_d[:, 0]))]
    unique_rows, indices = np.unique(t_d[:, :2], axis=0, return_inverse=True)
    grouped_values = np.split(t_d[:, 2], np.cumsum(np.bincount(indices))[:-1])
    grouped_values = [np.sort(np.unique(a)) for a in grouped_values]

    lines_s = np.concatenate(
        [
            hexLines_direction(
                a=diff_y * 2 / 3,
                i=diff_x,
                beta=beta - 0.1,
                off=unique_rows[i],
                directions=grouped_values[i].astype(int),
            )
            for i in range(len(unique_rows))
        ]
    )
    for l in lines_s:
        plt.plot(*l.transpose(), "w-", lw=lw, scalex=False, scaley=False, color="white")

    lines_s = np.concatenate(
        [
            hexLines_direction(
                a=diff_y * 2 / 3,
                i=diff_x,
                beta=beta,
                off=unique_rows[i],
                directions=grouped_values[i].astype(int),
            )
            for i in range(len(unique_rows))
        ]
    )
    for l in lines_s:
        plt.plot(*l.transpose(), "w-", lw=lw, scalex=False, scaley=False, color=color)


def hexLines(a=None, i=None, off=[0, 0]):
    """regular hexagon segment lines as `(xy1,xy2)` in clockwise
    order with points in line sorted top to bottom
    for irregular hexagon pass both `a` (vertical) and `i` (horizontal)"""
    if a is None:
        a = 2 / np.sqrt(3) * i
    if i is None:
        i = np.sqrt(3) / 2 * a
    h = a / 2
    xy = np.array(
        [
            [[0, a], [i, h]],
            [[i, h], [i, -h]],
            [[i, -h], [0, -a]],
            [[-i, -h], [0, -a]],  # flipped
            [[-i, h], [-i, -h]],  # flipped
            [[0, a], [-i, h]],  # flipped
        ]
    )
    return xy + off


def hexLines_direction(a=None, i=None, beta=0.88, off=[0, 0], directions=[1, 2, 5]):
    """regular hexagon segment lines as `(xy1,xy2)` in clockwise
    order with points in line sorted top to bottom
    for irregular hexagon pass both `a` (vertical) and `i` (horizontal)"""
    if a is None:
        a = 2 / np.sqrt(3) * i
    if i is None:
        i = np.sqrt(3) / 2 * a
    h = a / 2
    xy = np.array(
        [
            [[0, a], [i, h]],
            [[i, h], [i, -h]],
            [[i, -h], [0, -a]],
            [[-i, -h], [0, -a]],  # flipped
            [[-i, h], [-i, -h]],  # flipped
            [[0, a], [-i, h]],  # flipped
        ]
    )
    xy_beta = xy * beta
    xy_extents = np.array(
        [get_extent_points(i, xy, xy_beta) for i in range(len(xy_beta))]
    )
    for i_di in range(len(directions)):
        i_di_left = i_di - 1
        i_di_right = (i_di + 1) % len(directions)
        if (directions[i_di_left] + 1) % 6 != directions[i_di]:
            if directions[i_di] < 3:
                xy_beta[directions[i_di]][0] = xy_extents[directions[i_di]][0]
            else:
                xy_beta[directions[i_di]][1] = xy_extents[directions[i_di]][0]
        if (directions[i_di_right] - 1) % 6 != directions[i_di]:
            if directions[i_di] < 3:
                xy_beta[directions[i_di]][1] = xy_extents[directions[i_di]][1]
            else:
                xy_beta[directions[i_di]][0] = xy_extents[directions[i_di]][1]
    xy_beta = xy_beta[directions]
    return xy_beta + off


def get_extent_points(i, xy, xy_beta):
    x1 = xy_beta[i][0][0]
    y1 = xy_beta[i][0][1]
    x2 = xy_beta[i][1][0]
    y2 = xy_beta[i][1][1]

    if i == 0:
        j = 5
    else:
        j = i - 1
    x3 = xy[j][0][0]
    y3 = xy[j][0][1]
    x4 = xy[j][1][0]
    y4 = xy[j][1][1]

    if i == 5:
        k = 0
    else:
        k = i + 1
    x5 = xy[k][0][0]
    y5 = xy[k][0][1]
    x6 = xy[k][1][0]
    y6 = xy[k][1][1]

    x_left, y_left = find_intersection(x1, y1, x2, y2, x3, y3, x4, y4)
    x_right, y_right = find_intersection(x1, y1, x2, y2, x5, y5, x6, y6)

    return [[x_left, y_left], [x_right, y_right]]


def find_intersection(x1, y1, x2, y2, x3, y3, x4, y4):
    if np.isinf((y2 - y1) / (x2 - x1)):
        # Line 1 is vertical
        x_intersect = x1
        slope2 = (y4 - y3) / (x4 - x3)
        intercept2 = y3 - slope2 * x3
        y_intersect = slope2 * x_intersect + intercept2
    elif np.isinf((y4 - y3) / (x4 - x3)):
        # Line 2 is vertical
        x_intersect = x3
        slope1 = (y2 - y1) / (x2 - x1)
        intercept1 = y1 - slope1 * x1
        y_intersect = slope1 * x_intersect + intercept1
    else:
        # Neither line is vertical, calculate intersection normally
        slope1 = (y2 - y1) / (x2 - x1)
        slope2 = (y4 - y3) / (x4 - x3)
        intercept1 = y1 - slope1 * x1
        intercept2 = y3 - slope2 * x3
        x_intersect = (intercept2 - intercept1) / (slope1 - slope2)
        y_intersect = slope1 * x_intersect + intercept1
    # Return the coordinates of the intersection point
    return x_intersect, y_intersect


def getDirections(off=[0, 0]):
    xy = np.array(
        [
            [off[0], off[1], 0],
            [off[0], off[1], 1],
            [off[0], off[1], 2],
            [off[0], off[1], 3],
            [off[0], off[1], 4],
            [off[0], off[1], 5],
        ]
    )
    return xy


def get_nearest_index(array1, array2):
    distances = cdist(array1, array2)  # Calculate pairwise distances
    nearest_indices_c2 = np.argmin(distances, axis=1)
    return nearest_indices_c2


def plot_outline_only_1(
    verts, pa, tpa, plt, diff_y, diff_x, lw, color="white", beta=0.88
):
    nv = verts[pa == tpa]
    lines = np.concatenate(
        [hexLines(a=diff_y * 2 / 3, i=diff_x, off=off) for off in nv]
    )
    uls, c = np.unique(lines.round(2), axis=0, return_counts=True)
    for l in uls[c == 1]:
        plt.plot(*l.transpose(), "w-", lw=lw, scalex=False, scaley=False, color=color)


def dim_plot_cluster(
    adata_path, annotation_table_path, color_table_path, out_folder_path
):
    adata = sc.read(adata_path)
    if not os.path.exists(out_folder_path):
        os.mkdir(out_folder_path)
    annotation_table = pd.read_csv(annotation_table_path, header=0, na_filter=False)
    annotation_table["Pathology_annotations"] = (
        annotation_table["Pathology_annotations"].str.strip().replace("", "Unlabelled")
    )
    tempObs = adata.obs.merge(
        annotation_table, left_index=True, right_on="Barcode", how="left"
    )
    tempObs.set_index("Barcode", inplace=True)
    adata.obs = tempObs
    unique_PAs = annotation_table["Pathology_annotations"].unique()
    mapping_array = pd.Series([i + 1 for i in range(len(unique_PAs))], index=unique_PAs)
    adata.obs["pa_index"] = adata.obs["Pathology_annotations"].map(mapping_array)
    color_df = pd.read_csv(color_table_path, sep="\t")
    spot_array = np.array(
        [
            [x, y]
            for x, y in zip(
                (max(adata.obs["y_array"]) - adata.obs["y_array"]),
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            )
        ]
    )
    plt.figure(figsize=(9, 7), linewidth=0.0)
    hb = plt.hexbin(
        max(adata.obs["y_array"]) - adata.obs["y_array"],
        2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
        C=adata.obs["refined_pred"],
        gridsize=(
            math.floor((max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2),
            math.floor((max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2),
        ),
        edgecolors="face",
        linewidth=0.1,
        cmap=pylab.cm.get_cmap("tab20", len(np.unique(adata.obs["refined_pred"]))),
    )
    plt.axis("off")
    plt.gca().set_frame_on(False)
    cb = plt.colorbar(hb, shrink=0.6)
    cb.set_ticks(
        [
            (i + 0.5)
            * (len(np.unique(adata.obs["refined_pred"])) - 1)
            / len(np.unique(adata.obs["refined_pred"]))
            for i in range(len(np.unique(adata.obs["refined_pred"])))
        ]
    )
    cb.set_ticklabels(np.unique(adata.obs["refined_pred"]))
    plt.tight_layout()
    figure_path = "{0}/{1}.pdf".format(out_folder_path, "dim_plot_cluster")
    plt.savefig(figure_path)
    plt.close()


def dim_plot_cluster_ha(
    adata_path, annotation_table_path, color_table_path, out_folder_path
):
    adata = sc.read(adata_path)
    if not os.path.exists(out_folder_path):
        os.mkdir(out_folder_path)
    annotation_table = pd.read_csv(annotation_table_path, header=0, na_filter=False)
    annotation_table["Pathology_annotations"] = (
        annotation_table["Pathology_annotations"].str.strip().replace("", "Unlabelled")
    )
    tempObs = adata.obs.merge(
        annotation_table, left_index=True, right_on="Barcode", how="left"
    )
    tempObs.set_index("Barcode", inplace=True)
    adata.obs = tempObs
    unique_PAs = annotation_table["Pathology_annotations"].unique()
    mapping_array = pd.Series([i + 1 for i in range(len(unique_PAs))], index=unique_PAs)
    adata.obs["pa_index"] = adata.obs["Pathology_annotations"].map(mapping_array)
    color_df = pd.read_csv(color_table_path, sep="\t")
    spot_array = np.array(
        [
            [x, y]
            for x, y in zip(
                (max(adata.obs["y_array"]) - adata.obs["y_array"]),
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            )
        ]
    )
    plt.figure(figsize=(9, 7), linewidth=0.0)
    hb = plt.hexbin(
        max(adata.obs["y_array"]) - adata.obs["y_array"],
        2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
        C=adata.obs["refined_pred"],
        gridsize=(
            math.floor((max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2),
            math.floor((max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2),
        ),
        edgecolors="face",
        linewidth=0.1,
        cmap=pylab.cm.get_cmap("tab20", len(np.unique(adata.obs["refined_pred"]))),
    )
    plt.axis("off")
    plt.gca().set_frame_on(False)
    cb = plt.colorbar(hb, shrink=0.6)
    cb.set_ticks(
        [
            (i + 0.5)
            * (len(np.unique(adata.obs["refined_pred"])) - 1)
            / len(np.unique(adata.obs["refined_pred"]))
            for i in range(len(np.unique(adata.obs["refined_pred"])))
        ]
    )
    cb.set_ticklabels(np.unique(adata.obs["refined_pred"]))
    plt.tight_layout()
    verts = hb.get_offsets()
    unique_y = np.unique(np.array([item[1] for item in verts]))
    diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
    unique_x = np.unique(np.array([item[0] for item in verts]))
    diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
    index = get_nearest_index(verts, spot_array)
    pa_indexed = adata.obs.index[index]
    index_df = pd.Index(pa_indexed)
    pa = (
        annotation_table.loc[annotation_table["Barcode"].isin(pa_indexed)]
        .set_index("Barcode")
        .reindex(index_df)["Pathology_annotations"]
        .str.strip()
    )
    for tpa in np.unique(pa):
        if tpa != "":
            tcolor = color_df.loc[color_df["Tissue"] == tpa, "Color"].values[0]
            plot_outline(
                verts=verts,
                pa=pa,
                tpa=tpa,
                plt=plt,
                diff_y=diff_y,
                diff_x=diff_x,
                lw=0.5,
                beta=0.88,
                color=tcolor,
            )
    figure_path = "{0}/{1}.pdf".format(out_folder_path, "dim_plot_cluster_ha")
    plt.savefig(figure_path)
    plt.close()


def dim_plot_ha(adata_path, annotation_table_path, color_table_path, out_folder_path):
    adata = sc.read(adata_path)
    if not os.path.exists(out_folder_path):
        os.mkdir(out_folder_path)
    annotation_table = pd.read_csv(annotation_table_path, header=0, na_filter=False)
    annotation_table["Pathology_annotations"] = (
        annotation_table["Pathology_annotations"].str.strip().replace("", "Unlabelled")
    )
    tempObs = adata.obs.merge(
        annotation_table, left_index=True, right_on="Barcode", how="left"
    )
    tempObs.set_index("Barcode", inplace=True)
    adata.obs = tempObs
    unique_PAs = annotation_table["Pathology_annotations"].unique()
    mapping_array = pd.Series([i + 1 for i in range(len(unique_PAs))], index=unique_PAs)
    adata.obs["pa_index"] = adata.obs["Pathology_annotations"].map(mapping_array)
    color_df = pd.read_csv(color_table_path, sep="\t")
    color_values = [
        color_df.loc[color_df["Tissue"] == t, "Color"].iloc[0] for t in unique_PAs
    ]
    color_cmap = mp.colors.ListedColormap(color_values)
    spot_array = np.array(
        [
            [x, y]
            for x, y in zip(
                (max(adata.obs["y_array"]) - adata.obs["y_array"]),
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            )
        ]
    )
    plt.figure(figsize=(9, 7), linewidth=0.0)
    hb = plt.hexbin(
        max(adata.obs["y_array"]) - adata.obs["y_array"],
        2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
        C=adata.obs["pa_index"],
        gridsize=(
            math.floor((max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2),
            math.floor((max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2),
        ),
        edgecolors="face",
        linewidth=0.1,
        cmap=color_cmap,
    )
    plt.axis("off")
    plt.gca().set_frame_on(False)
    cb = plt.colorbar(hb, shrink=0.6)
    cb.set_ticks(
        [
            (i + 1.5) * (len(unique_PAs) - 1) / len(unique_PAs)
            for i in range(len(unique_PAs))
        ]
    )
    cb.set_ticklabels(unique_PAs)
    plt.tight_layout()
    figure_path = "{0}/{1}.pdf".format(out_folder_path, "dim_plot_ha")
    plt.savefig(figure_path)
    plt.close()


def dim_plot_ha_cluster(
    adata_path, annotation_table_path, color_table_path, out_folder_path
):
    adata = sc.read(adata_path)
    if not os.path.exists(out_folder_path):
        os.mkdir(out_folder_path)
    annotation_table = pd.read_csv(annotation_table_path, header=0, na_filter=False)
    annotation_table["Pathology_annotations"] = (
        annotation_table["Pathology_annotations"].str.strip().replace("", "Unlabelled")
    )
    tempObs = adata.obs.merge(
        annotation_table, left_index=True, right_on="Barcode", how="left"
    )
    tempObs.set_index("Barcode", inplace=True)
    adata.obs = tempObs
    unique_PAs = annotation_table["Pathology_annotations"].unique()
    mapping_array = pd.Series([i + 1 for i in range(len(unique_PAs))], index=unique_PAs)
    adata.obs["pa_index"] = adata.obs["Pathology_annotations"].map(mapping_array)
    color_df = pd.read_csv(color_table_path, sep="\t")
    color_values = [
        color_df.loc[color_df["Tissue"] == t, "Color"].iloc[0] for t in unique_PAs
    ]
    color_cmap = mp.colors.ListedColormap(color_values)
    spot_array = np.array(
        [
            [x, y]
            for x, y in zip(
                (max(adata.obs["y_array"]) - adata.obs["y_array"]),
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            )
        ]
    )
    plt.figure(figsize=(9, 7), linewidth=0.0)
    hb = plt.hexbin(
        max(adata.obs["y_array"]) - adata.obs["y_array"],
        2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
        C=adata.obs["pa_index"],
        gridsize=(
            math.floor((max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2),
            math.floor((max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2),
        ),
        edgecolors="face",
        linewidth=0.1,
        cmap=color_cmap,
    )
    plt.axis("off")
    plt.gca().set_frame_on(False)
    cb = plt.colorbar(hb, shrink=0.6)
    cb.set_ticks(
        [
            (i + 1.5) * (len(unique_PAs) - 1) / len(unique_PAs)
            for i in range(len(unique_PAs))
        ]
    )
    cb.set_ticklabels(unique_PAs)
    plt.tight_layout()
    verts = hb.get_offsets()
    unique_y = np.unique(np.array([item[1] for item in verts]))
    diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
    unique_x = np.unique(np.array([item[0] for item in verts]))
    diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
    index = get_nearest_index(verts, spot_array)
    cluster_indexed = adata.obs.index[index]
    cluster_index_df = adata.obs.loc[cluster_indexed, "refined_pred"]

    cmap_cluster = pylab.cm.get_cmap("tab20", len(np.unique(adata.obs["refined_pred"])))
    hex_colors_cluster = [
        mp.colors.rgb2hex(color)
        for color in cmap_cluster(
            [i for i in range(len(np.unique(adata.obs["refined_pred"])))]
        )
    ]
    for tci in np.unique(cluster_index_df):
        if tci != "":
            tcolor = hex_colors_cluster[tci]
            plot_outline(
                verts=verts,
                pa=cluster_index_df,
                tpa=tci,
                plt=plt,
                diff_y=diff_y,
                diff_x=diff_x,
                lw=0.5,
                beta=0.88,
                color=tcolor,
            )
    figure_path = "{0}/{1}.pdf".format(out_folder_path, "dim_plot_ha_cluster")
    plt.savefig(figure_path)
    plt.close()


def feature_plot_hex_cluster(
    adata_path, marker_table_path, out_folder_path, marker_index=0
):
    if isinstance(adata_path, str):
        adata = sc.read(adata_path)
    else:
        adata = adata_path

    if not os.path.exists(out_folder_path):
        os.mkdir(out_folder_path)
    marker_table = pd.read_csv(marker_table_path, header=0, na_filter=False, sep="\t")
    marker_table.columns.values[marker_index] = "Marker"
    spot_array = np.array(
        [
            [x, y]
            for x, y in zip(
                (max(adata.obs["y_array"]) - adata.obs["y_array"]),
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            )
        ]
    )
    for marker in marker_table["Marker"].unique():
        if marker in adata.var_names:
            plt.figure(figsize=(8, 7), linewidth=0.0)
            hb = plt.hexbin(
                max(adata.obs["y_array"]) - adata.obs["y_array"],
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
                C=adata[adata.obs.index, marker].X.toarray().flatten(),
                gridsize=(
                    math.floor(
                        (max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2
                    ),
                    math.floor(
                        (max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2
                    ),
                ),
                edgecolors="face",
                linewidth=0.1,
                cmap="gnuplot2",
            )
            plt.axis("off")
            plt.gca().set_frame_on(False)
            cb = plt.colorbar(hb, shrink=0.6)
            cb.set_label("log1pEXP")
            plt.title(marker)
            plt.tight_layout()
            verts = hb.get_offsets()
            unique_y = np.unique(np.array([item[1] for item in verts]))
            diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
            unique_x = np.unique(np.array([item[0] for item in verts]))
            diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
            index = get_nearest_index(verts, spot_array)
            cluster_indexed = adata.obs.index[index]
            cluster_index_df = adata.obs.loc[cluster_indexed, "refined_pred"]
            cmap_cluster = pylab.cm.get_cmap(
                "tab20", len(np.unique(adata.obs["refined_pred"]))
            )
            hex_colors_cluster = [
                mp.colors.rgb2hex(color)
                for color in cmap_cluster(
                    [i for i in range(len(np.unique(adata.obs["refined_pred"])))]
                )
            ]
            for tci in np.unique(cluster_index_df):
                if tci != "":
                    tcolor = hex_colors_cluster[tci]
                    plot_outline(
                        verts=verts,
                        pa=cluster_index_df,
                        tpa=tci,
                        plt=plt,
                        diff_y=diff_y,
                        diff_x=diff_x,
                        lw=0.5,
                        beta=0.88,
                        color=tcolor,
                    )
            figure_path = "{0}/{1}_cluster.pdf".format(out_folder_path, marker)
            plt.savefig(figure_path)
            plt.close()


def get_module_score(
    adata_path,
    marker_table_path,
    annotation_table_path,
    color_table_path,
    out_folder_path,
):
    if isinstance(adata_path, str):
        adata = sc.read(adata_path)
    else:
        adata = adata_path

    # adata = sc.read(
    #     "/rsrch3/scratch/genomic_med/ychu2/projects/project26/result/PM3131/3_SpaGCN/HKU01/alpha_1.5_beta_30/p_0.4_nClusters_12/domain.h5ad"
    # )

    annotation_table = pd.read_csv(annotation_table_path, header=0, na_filter=False)
    annotation_table["Pathology_annotations"] = (
        annotation_table["Pathology_annotations"].str.strip().replace("", "Unlabelled")
    )
    tempObs = adata.obs.merge(
        annotation_table, left_index=True, right_on="Barcode", how="left"
    )
    tempObs.set_index("Barcode", inplace=True)
    adata.obs = tempObs
    unique_PAs = annotation_table["Pathology_annotations"].unique()
    mapping_array = pd.Series([i + 1 for i in range(len(unique_PAs))], index=unique_PAs)
    adata.obs["pa_index"] = adata.obs["Pathology_annotations"].map(mapping_array)
    color_df = pd.read_csv(color_table_path, sep="\t")
    spot_array = np.array(
        [
            [x, y]
            for x, y in zip(
                (max(adata.obs["y_array"]) - adata.obs["y_array"]),
                2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            )
        ]
    )

    marker_table = pd.read_csv(marker_table_path, header=0, na_filter=False, sep="\t")
    # marker_table = pd.read_csv(
    #     "/rsrch3/scratch/genomic_med/ychu2/projects/project26/knowledge/private/MPs41-ITH.txt",
    #     header=0,
    #     na_filter=False,
    #     sep="\t",
    # )

    for cn in marker_table.columns.tolist():
        gene_list = marker_table[cn].tolist()
        sc.tl.score_genes(adata, gene_list, score_name=cn)

        plt.figure(figsize=(8, 7), linewidth=0.0)
        hb = plt.hexbin(
            max(adata.obs["y_array"]) - adata.obs["y_array"],
            2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            C=adata.obs[cn],
            gridsize=(
                math.floor((max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2),
                math.floor((max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2),
            ),
            edgecolors="face",
            linewidth=0.1,
            cmap="gnuplot2",
        )
        plt.axis("off")
        plt.gca().set_frame_on(False)
        cb = plt.colorbar(hb, shrink=0.6)
        cb.set_label("Avg log1pEXP")
        plt.title(cn)
        plt.tight_layout()
        verts = hb.get_offsets()
        unique_y = np.unique(np.array([item[1] for item in verts]))
        diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
        unique_x = np.unique(np.array([item[0] for item in verts]))
        diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
        index = get_nearest_index(verts, spot_array)
        pa_indexed = adata.obs.index[index]
        index_df = pd.Index(pa_indexed)
        pa = (
            annotation_table.loc[annotation_table["Barcode"].isin(pa_indexed)]
            .set_index("Barcode")
            .reindex(index_df)["Pathology_annotations"]
            .str.strip()
        )
        for tpa in np.unique(pa):
            if tpa != "":
                tcolor = color_df.loc[color_df["Tissue"] == tpa, "Color"].values[0]
                plot_outline(
                    verts=verts,
                    pa=pa,
                    tpa=tpa,
                    plt=plt,
                    diff_y=diff_y,
                    diff_x=diff_x,
                    lw=0.5,
                    beta=0.88,
                    color=tcolor,
                )
        figure_path = "{0}/{1}.pdf".format(out_folder_path, cn)
        plt.savefig(figure_path)
        plt.close()

        plt.figure(figsize=(8, 7), linewidth=0.0)
        hb = plt.hexbin(
            max(adata.obs["y_array"]) - adata.obs["y_array"],
            2 * (max(adata.obs["x_array"]) - adata.obs["x_array"]),
            C=adata.obs[cn],
            gridsize=(
                math.floor((max(adata.obs["y_array"]) - min(adata.obs["y_array"])) / 2),
                math.floor((max(adata.obs["x_array"]) - min(adata.obs["x_array"])) / 2),
            ),
            edgecolors="face",
            linewidth=0.1,
            cmap="gnuplot2",
        )
        plt.axis("off")
        plt.gca().set_frame_on(False)
        cb = plt.colorbar(hb, shrink=0.6)
        cb.set_label("log1pEXP")
        plt.title(cn)
        plt.tight_layout()
        verts = hb.get_offsets()
        unique_y = np.unique(np.array([item[1] for item in verts]))
        diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
        unique_x = np.unique(np.array([item[0] for item in verts]))
        diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
        index = get_nearest_index(verts, spot_array)
        cluster_indexed = adata.obs.index[index]
        cluster_index_df = adata.obs.loc[cluster_indexed, "refined_pred"]
        cmap_cluster = pylab.cm.get_cmap(
            "tab20", len(np.unique(adata.obs["refined_pred"]))
        )
        hex_colors_cluster = [
            mp.colors.rgb2hex(color)
            for color in cmap_cluster(
                [i for i in range(len(np.unique(adata.obs["refined_pred"])))]
            )
        ]
        for tci in np.unique(cluster_index_df):
            if tci != "":
                tcolor = hex_colors_cluster[tci]
                plot_outline(
                    verts=verts,
                    pa=cluster_index_df,
                    tpa=tci,
                    plt=plt,
                    diff_y=diff_y,
                    diff_x=diff_x,
                    lw=0.5,
                    beta=0.88,
                    color=tcolor,
                )
        figure_path = "{0}/{1}_cluster.pdf".format(out_folder_path, cn)
        plt.savefig(figure_path)
        plt.close()

    return adata


def main():
    import argparse

    parser = argparse.ArgumentParser(description="plot gene expression")

    parser.add_argument("--adataPath", dest="adataPath", help="adata path")
    parser.add_argument(
        "--markerTablePath", dest="markerTablePath", help="marker table path"
    )
    parser.add_argument(
        "--annotationTablePath",
        dest="annotationTablePath",
        help="annotation table path",
    )
    parser.add_argument(
        "--colorTablePath", dest="colorTablePath", help="color table path"
    )
    parser.add_argument("--outFolderPath", dest="outFolderPath", help="out folder path")
    parser.add_argument(
        "--markerIndex", dest="markerIndex", type=int, default=0, help="marker index"
    )

    args = parser.parse_args()

    print(args)

    # dim_plot_cluster(adata_path=args.adataPath,
    #                  annotation_table_path=args.annotationTablePath,
    #                  color_table_path=args.colorTablePath,
    #                  out_folder_path=os.path.dirname(args.outFolderPath))

    # dim_plot_ha(adata_path=args.adataPath,
    #             annotation_table_path=args.annotationTablePath,
    #             color_table_path=args.colorTablePath,
    #             out_folder_path=os.path.dirname(args.outFolderPath))

    # dim_plot_cluster_ha(adata_path=args.adataPath,
    #                     annotation_table_path=args.annotationTablePath,
    #                     color_table_path=args.colorTablePath,
    #                     out_folder_path=os.path.dirname(args.outFolderPath))

    # dim_plot_ha_cluster(adata_path=args.adataPath,
    #                     annotation_table_path=args.annotationTablePath,
    #                     color_table_path=args.colorTablePath,
    #                     out_folder_path=os.path.dirname(args.outFolderPath))

    # feature_plot_hex_ha(
    #     adata_path=args.adataPath,
    #     marker_table_path=args.markerTablePath,
    #     annotation_table_path=args.annotationTablePath,
    #     color_table_path=args.colorTablePath,
    #     out_folder_path=args.outFolderPath,
    #     marker_index=args.markerIndex,
    # )

    # feature_plot_hex_cluster(
    #     adata_path=args.adataPath,
    #     marker_table_path=args.markerTablePath,
    #     out_folder_path=args.outFolderPath,
    #     marker_index=args.markerIndex,
    # )

    get_module_score(
        adata_path=args.adataPath,
        marker_table_path=args.markerTablePath,
        out_folder_path=args.outFolderPath,
        annotation_table_path=args.annotationTablePath,
        color_table_path=args.colorTablePath,
    )


if __name__ == "__main__":
    main()
