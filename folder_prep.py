#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 13:26:44 2022

@author: jchu
"""

import scanpy as sc
import numpy as np
from scanpy import read_10x_h5
import os,csv,re
import shutil
import pandas as pd 
# local file
from paste_utils import prepare_paste_input_dir, load_slices, plot_multiple_expression
import paste as pst
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from copy import deepcopy
from anndata import AnnData
from typing import Optional
from matplotlib import image
import matplotlib.image as mpimg
import json
import pickle
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#inatll opencv in python
#!pip install opencv-python
import cv2


def prefilter_genes(adata,min_counts=None,max_counts=None,min_cells=10,max_cells=None):
    if min_cells is None and min_counts is None and max_cells is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp=np.asarray([True]*adata.shape[1],dtype=bool)
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_cells=min_cells)[0]) if min_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_cells=max_cells)[0]) if max_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_counts=min_counts)[0]) if min_counts is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_counts=max_counts)[0]) if max_counts is not None  else id_tmp
    adata._inplace_subset_var(id_tmp)

def plot_cluster_per_slice(rjdata_list, w, h, 
                           res, slice_names, donor_names,
                           n_ct,out_dir,
                           color_key="pred",suffix="",
                           plot_color=["#F56867","#FEB915","#C798EE","#59BE86",
                                        "#7495D3","#D1D1D1","#6D1A9C","#15821E",
                                        "#3A84E6","#997273","#787878","#DB4C6C",
                                        "#9E7A7A","#554236","#AF5F3C","#93796C",
                                        "#F9BD3F","#DAB370","#877F6C","#268785"],
                          save_assignment=False,prefix="ten_slices",dpi=200,figsize=(22,22)):
    cm = dict(zip([str(x) for x in range(n_ct)], plot_color[:n_ct]))
    fig, axs = plt.subplots(w, h, figsize=(22,22))
    i = 0
    for ax_i in axs.reshape(-1):
        #for i, RJDatai in enumerate(rjdata_list):
        if i >= len(rjdata_list):
            break
        RJDatai = rjdata_list[i]
        ct_i = RJDatai.obs[color_key].unique().tolist()
        ct_i.sort(key = int)
        RJDatai.uns[color_key+"_colors"] = [cm[x] for x in ct_i]
        ax  = sc.pl.scatter(RJDatai,
                            alpha=1,
                            x="pixel_y_align",
                            y="pixel_x_align",
                            color=color_key,
                            #color_map=cm,
                            show=False,size=100000/RJDatai.shape[0],
                            ax = ax_i)    
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        ax.axes.set_xlim([-1000, 1000])
        ax.axes.set_ylim([-1000, 1000])
        ax.set(title=donor_names[i]+"_"+slice_names[i])
    #     axs[x,y] = ax
    #     plt.savefig("{dir}/segementation_joint_{res}_{tissueID_ref}.png".format(dir = out_dir, res = res, tissueID_ref = slice_names[i]), dpi=600)
    #     plt.close()
        if save_assignment:
            RJDatai.obs.to_csv("{dir}/obsdata_joint_{res}_{tissueID_ref}.csv".format(dir = out_dir,res = res, tissueID_ref =slice_names[i]))
        i += 1 

    plt.savefig("{dir}/{prefix}_spagcn_{res}{suffix}.png".format(dir = out_dir, prefix=prefix, res = res,suffix=suffix), dpi=dpi)
    plt.close()
    

#Folder preparation

# prepare_paste_input_dir("/mnt/depts/dept04/compbio/users/ezhao/data-raw/",
#                         "../animals/animal45/", ["106A","107B","108A"])

# prepare_paste_input_dir("/mnt/depts/dept04/compbio/users/ezhao/data-raw/", 
#                         "../animals/animal43/", ["072D","072A","074D"])

# prepare_paste_input_dir("/mnt/depts/dept04/compbio/users/ezhao/data-raw/", 
#                         "../animals/animal46/", ["075D","076C","074C"])    

prepare_paste_input_dir("/mnt/depts/dept04/compbio/users/ezhao/data-raw/",
                        "../animals/animal14/", ["075B","076D","072B"])

prepare_paste_input_dir("/mnt/depts/dept04/compbio/users/ezhao/data-raw/", 
                        "../animals/animal53/", ["106C","107C","108D","076A"])

prepare_paste_input_dir("/mnt/depts/dept04/compbio/users/ezhao/data-raw/", 
                        "../animals/animal13/", ["106B","107A","108C"])    
