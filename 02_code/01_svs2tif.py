# -*- coding: utf-8 -*-
# __author__ = "Pingjun Chen"
# __email__ =  "chenpingjun@gmx.com"
module load python/3.7.3-anaconda
conda create -n openslide_env
conda activate openslide_env
conda install -c conda-forge openslide
python
import openslide

import os, sys, pdb
import numpy as np
from scipy import misc
import matplotlib.pyplot as plt
from skimage import io
import openslide
import argparse


def save_svs_img(slide_filename, tile_size=8192):
	slide_file = openslide.OpenSlide(slide_filename)
	slide_width, slide_height = slide_file.dimensions

	# tile_arr = []
	slide_img = np.zeros((slide_height, slide_width, 3), np.uint8)	
	x_tile_num = int(np.floor((slide_width-1)/tile_size)) + 1
	y_tile_num = int(np.floor((slide_height-1)/tile_size)) + 1
	for iy in range(y_tile_num):	
		for ix in range(x_tile_num):
			start_x = ix * tile_size
			len_x = tile_size if (ix + 1) * tile_size < slide_width else (slide_width - start_x) 
			start_y = iy * tile_size
			len_y = tile_size if (iy + 1) * tile_size < slide_height else (slide_height - start_y)
			# tile_arr.append(((start_x, start_y), (len_x, len_y)))
			cur_tile = slide_file.read_region(location=(start_x, start_y), level=0, size=(len_x, len_y))
			slide_img[start_y:start_y+len_y, start_x:start_x+len_x, :] = np.array(cur_tile)[:,:,:3]

	slide_savename = os.path.splitext(slide_filename)[0] + '.tif'
	# misc.imsave(slide_savename, slide_img)
	io.imsave(slide_savename, slide_img)


def batch_convert_svs(slide_dir, tile_size=8192):
	svs_list = [svs_file for svs_file in os.listdir(slide_dir) if svs_file.endswith(".svs")]
	for ind, svs_file in enumerate(svs_list):
		print("Processing {}/{}: {}".format(ind+1, len(svs_list), svs_file))
		save_svs_img(svs_file, tile_size)


def set_args():
    parser = argparse.ArgumentParser(description="load svs")

    parser.add_argument("--filedir",        type=str, default="./")
    parser.add_argument("--filename",       type=str, default="124180.svs")
    parser.add_argument("--tile-size",      type=int, default=4096)    

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = set_args()

    # save_svs_img(args.filename, args.tile_size)
    batch_convert_svs(args.filedir, args.tile_size)