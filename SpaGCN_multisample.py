#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:15:13 2022

@author: stjhc
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

# prepare_paste_input_dir("./", "./example/", ["106C","106A","072A"])

# prepare_paste_input_dir("../spacerange/", "../animals/animal53/", ["106C","107C","108D","076A"])

# prepare_paste_input_dir("../spacerange/", "../animals/animal13/", ["106B","107A","108C"])    


#Center slice identification

# Load Slices
# Control
# data_dirs = ["./animals/animal14/"]*3 + ["./animals/animal53/"]*4 + ["./animals/animal13/"]*3
# slice_names = ["072B", "075B","076D"] + ["106C","107C","108D","076A"] + ["106B","107A","108C"]

# Cuprizone
# data_dirs = ["./animals/animal45/"]*3 + ["./animals/animal43/"]*3 + ["./animals/animal46/"]*3
# slice_names = ["106A","107B","108A"] + ["072D","072A","074D"] + ["075D","076C","074C"]

# All samples
data_dirs = ["./animals/animal14/"]*3 + ["./animals/animal53/"]*4 + ["./animals/animal13/"]*3 + ["./animals/animal45/"]*3 + ["./animals/animal43/"]*3 + ["./animals/animal46/"]*3+ ["./animals/animal39/"]*3 + ["./animals/animal40/"]*3 + ["./animals/animal41/"]*3
            
slice_names = ["072B", "075B","076D"] + ["106C","107C","108D","076A"] + ["106B","107A","108C"] + ["106A","107B","108A"] + ["072D","072A","074D"] + ["075D","076C","074C"] + ["106D","107D","108B"] + ["074B","072C","074A"] + ["075C","076B","075A"]


slices = load_slices(data_dirs,slice_names)

# Construct a center slice
## choose one of the slices as the coordinate reference for the center slice,
## i.e. the center slice will have the same number of spots as this slice and
## the same coordinates.
slice1 = slices[0]
initial_slice = slice1.copy()
lmbda = len(slices)*[1/len(slices)]

## Possible to pass in an initial pi (as keyword argument pis_init) 
## to improve performance, see Tutorial.ipynb notebook for more details.pooooooooooooooo-0
center_slice, pis = pst.center_align(initial_slice, slices, lmbda) 

## The low dimensional representation of our center slice is held 
## in the matrices W and H, which can be used for downstream analyses
W = center_slice.uns['paste_W']
H = center_slice.uns['paste_H']

#Alignment to center slice

center, new_slices_all = pst.stack_slices_center(center_slice, slices, pis)

center_color = 'orange'
#slices_colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3']
palette = sns.color_palette(None, len(slices))
palette

plt.figure(figsize=(7,7))
#pst.plot_slice(center,center_color,s=10)
for i in range(len(new_slices_all)):
    slice_i = deepcopy(new_slices_all[i])
    # reversing x-y coordinate
    slice_i.obsm['spatial'] = slice_i.obsm['spatial'][:, [1, 0]]
    pst.plot_slice(slice_i,palette[i],s=10)

handle_list = [mpatches.Patch(color=palette[i], label=slice_names[i]) for i in range(len(palette))]
plt.legend(handles=handle_list)
plt.gca().invert_yaxis()
plt.axis('off')
plt.show()

# Save aligned result
out_dir = "./paste/integration_test/"
os.makedirs(out_dir, exist_ok=True)
#file_path = "{out}/ten_slices_alignment_res.pl".format(out=out_dir)
file_path = "{out}/all_slices_alignment_res.pl".format(out=out_dir)
#file_path = "{out}/cup_slices_alignment_res.pl".format(out=out_dir)

# data_dirs = ["./animals/animal14/"]*3 + ["./animals/animal53/"]*4 + ["./animals/animal13/"]*3
# slice_names = ["072B", "075B","076D"] + ["106C","107C","108D","076A"] + ["106B","107A","108C"]

# data_dirs = ["../animals/animal14/"]*3 + ["../animals/animal53/"]*4 + ["../animals/animal13/"]*3
# slice_names = ["072B", "075B","076D"] + ["106C","107C","108D","076A"] + ["106B","107A","108C"]

with open(file_path, 'wb') as fp:
    pickle.dump([center, new_slices_all], fp)
    
#SpaGCN
# read saved result
in_dir = "./paste/integration_test/"
#file_path = "{out}/all_slices_alignment_res.pl".format(out=in_dir)

with open(file_path, 'rb') as fp:
    center, new_slices_all = pickle.load(fp)
    
# data_dirs = ["../animals/animal14/"]*3 + ["../animals/animal53/"]*4 + ["../animals/animal13/"]*3
# slice_names = ["072B", "075B","076D"] + ["106C","107C","108D","076A"] + ["106B","107A","108C"]
donor_names = ["Animal14"]*3+["Animal53"]*4+["Animal13"]*3+["Animal45"]*3+["Animal43"]*3+["Animal46"]*3+["Animal39"]*3+["Animal40"]*3+["Animal41"]*3
# donor_names = ["Animal45"]*3+["Animal43"]*3+["Animal46"]*3

#%timeit
suffix = "_p0d80_scaled"
out_dir = "./paste/integration_test{}/".format(suffix)
os.makedirs(out_dir, exist_ok=True)
nclusters = 17 # CAREFUL: the plotting function below only has 20 colors, remember to pass in more colors if ncluster >20
p=0.8  #p: percentage of total expression contributed by neighborhoods.
os.makedirs(out_dir, exist_ok=True)
#normalize gene expression 
slice_list = deepcopy(new_slices_all)

#ReferenceData_list = spatial_list = ReferenceData_list = [None]*len(slice_list)
adata_list = [None]*len(slice_list)
l_list = [None]*len(slice_list)
adj_list = [None]*len(slice_list)

for i in range(len(slice_list)):
    adata = deepcopy(slice_list[i])
    spg.prefilter_genes(adata,min_cells=20) # avoiding all gene is zeros
    adata.var_names_make_unique()
    adata.raw=adata
    sc.pp.normalize_per_cell(adata, min_counts=0)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    adata.obs['batch'] = slice_names[i]
    adata.obs['pixel_x_align'] = adata.obsm['spatial'][:,0].tolist()
    adata.obs['pixel_y_align'] = adata.obsm['spatial'][:,1].tolist()
    x = adata.obs['pixel_x_align']
    y = adata.obs['pixel_y_align']

    ### ------------------------------------------------------------------------------------------------------- ###
    ###        Step 2: Runing SpaGCN to obtain the spatial domains
    ### ------------------------------------------------------------------------------------------------------- ###
    adj=spg.calculate_adj_matrix(x=x,y=y,histology=False)
    adj_list[i] = adj
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    l_list[i] = l
    adata_list[i] = adata

#Set seed
r_seed=t_seed=n_seed=2022

def find_resolution_multispagcn(n_clusters): 
    #adata = adata_.copy()
    obtained_clusters = -1
    iteration = 0
    resolutions = [0., 1.]
    
    while obtained_clusters != n_clusters and iteration < 50:
        current_res = sum(resolutions)/2
        clf=spg.multiSpaGCN()
        # setting seeds! Need to explicitly set it for each source of randomness. 
        random.seed(r_seed)
        torch.manual_seed(t_seed)
        np.random.seed(n_seed)

        clf.train(deepcopy(adata_list),adj_list,l_list,init_spa=True,init="louvain",res=current_res, tol=5e-3, lr=0.05, max_epochs=200)
        y_pred, prob = clf.predict()
        labels = y_pred.astype('str')
        obtained_clusters = len(np.unique(labels))
        print(str(obtained_clusters))
        if obtained_clusters < n_clusters:
            resolutions[0] = current_res
        else:
            resolutions[1] = current_res
        
        iteration = iteration + 1
        print("res = {} found {} clusters".format(str(current_res),str(obtained_clusters)))
        
    return clf, current_res


#Multi-section SpaGCN
clf,res = find_resolution_multispagcn(nclusters)
#Get results
y_pred, prob = clf .predict()
adata_all = clf.adata_all
adata_all.obs["pred"] = y_pred
adata_all.obs["pred"] = adata_all.obs["pred"].astype('str')
adata_all.obs["pred"].value_counts().sort_index()

rjdata_list = [None]*len(slice_names)
for i in range(len(slice_names)):
    slice_name = slice_names[i]
    ref_id = adata_all.obs["batch"] == slice_name
    RJData=adata_all[ref_id,:].copy()
    #print(RJData.obs["pred"].value_counts().sort_index())
    rjdata_list[i] = RJData

# create subplots
w = 3
h = 3
n_ct = len(adata_all.obs["pred"].unique())

plot_cluster_per_slice(rjdata_list,w,h,res,slice_names,donor_names,n_ct,out_dir,suffix="_s2022{}".format(suffix))

df = adata_all.obs
df.to_csv("all_results_17clusters.csv", index = True, header=True)
