#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:17:09 2022

@author: jchu
"""

# ml load anaconda3
# conda activate st_env
# cd /home/jchu/st/spaGCN 

import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc, anndata as ad
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#inatll opencv in python
#!pip3 install opencv-python
import cv2
import json

from scanpy import read_10x_h5

    
block=[
       # '072A', '072B', '072C', '072D', 
       #  '074A', '074B', '074C', '074D',
       #  '075A', '075B', '075C', '075D',
       #  '076A', '076B', '076C', '076D',
       #  '106A', '106B', '106C', 
        '106D',
       '107A', '107B', '107C', '107D',
       '108A', '108B', '108C', '108D']
n_blocks=len(block)

#Current path is /home/jchu/st/spaGCN on edge
for i in block:
    adata=sc.read("./out/cuprizone/" + i + "_sample_data.h5ad")
    img=cv2.imread("/edgehpc/dept/compbio/users/ezhao/data-raw/" + i + "/spatial/tissue_hires_image.png")
 
    x_array=adata.obs["x_array"].tolist()
    y_array=adata.obs["y_array"].tolist()
    x_pixel=adata.obs["x_pixel"].tolist()
    y_pixel=adata.obs["y_pixel"].tolist()   

    f=open("/edgehpc/dept/compbio/users/ezhao/data-raw/" + i + "/spatial/scalefactors_json.json")
    scalefactor=json.load(f)
    #Adjust the image size by the scale factor
    x_pixel=[round(i * scalefactor['tissue_hires_scalef']) for i in x_pixel]
    y_pixel=[round(i * scalefactor['tissue_hires_scalef']) for i in y_pixel]
           

    img_new=img.copy()
    for j in range(len(x_pixel)):
        x=x_pixel[j]
        y=y_pixel[j]
        img_new[int(x-10):int(x+10), int(y-10):int(y+10),:]=255
    
    cv2.imwrite("./out/cuprizone/" + i + "_map.jpg", img_new)   
    
    
    adata.var_names_make_unique()
    spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    
    
    s=1
    b=49
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
    #calculate the adjacent matrix  
    
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    
    p=0.5 
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    
    #We set the number of clusters as 16 to compare with BayesSpace
    n_clusters=16
    #Set seed
    r_seed=t_seed=n_seed=100
    #Search for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    #Do cluster refinement(optional)
    #shape="hexagon" for Visium data, "square" for ST data.
    adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    #Save results
    #adata.write_h5ad("./out/cuprizone/" + i + "_results.h5ad")
    
    df = adata.obs
    df.to_csv("./out/cuprizone/" + i + "_results.csv", index = False, header=True)
    
    # adata=sc.read("./out/cuprizone/" + i + "_results.h5ad")
    # #Set colors used
    # plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
    # #Plot spatial domains
    # domains="pred"
    # num_celltype=len(adata.obs[domains].unique())
    # adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
    # ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
    # ax.set_aspect('equal', 'box')
    # ax.axes.invert_yaxis()
    # plt.savefig("./out/cuprizone/" + i + "_pred.png", dpi=600)
    # plt.close()
    
    # #Plot refined spatial domains
    # domains="refined_pred"
    # num_celltype=len(adata.obs[domains].unique())
    # adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
    # ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
    # ax.set_aspect('equal', 'box')
    # ax.axes.invert_yaxis()
    # plt.savefig("./out/cuprizone/" + i + "_refined_pred.png", dpi=600)
    # plt.close()