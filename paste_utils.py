import scanpy as sc
import numpy as np
from scanpy import read_10x_h5
import os 
import shutil
import pandas as pd 

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import paste as pst
import seaborn as sns
from typing import Optional
from matplotlib import image
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import json
import pickle
from copy import deepcopy
from anndata import AnnData

"""
    Organize files in formats accepted by paste algorithms
        1. Read filtered_feature_bc_matrix.h5
        2. Make gene names unique 
        3. load position list 
        4. keep position only for the relevant cells
        5. extract pixel.x and pixel.y column, save as column 1 and 2 in _coor.csv file 
        6. change variable name to be defined by the spot coordinate 
    Output:
        1. file_list[i].csv: storing spot-gene matrix 
        2. file_list[i]_coor.csv: x and y coordinate 

    Example input: prepare_paste_input_dir("../spacerange/", "../animals/animal14/", ["072B","075B","076D"])

"""

def prepare_paste_input_dir(in_dir, out_dir, file_list):

    for f in file_list:
        print("processing {}".format(f))
        adata = read_10x_h5("{}{}/filtered_feature_bc_matrix.h5".format(in_dir,f))
        adata.var_names_make_unique()
        pos = pd.read_csv("{}{}/spatial/tissue_positions_list.csv".format(in_dir,f),header=None)
        # select for barcodes present in filtered bc matrix
        common_bc = list(set(pos[0]).intersection(list(adata.obs_names)))
        # subset pos for the barcodes 
        pos.index = list(pos[0])
        pos_sel = pos.loc[common_bc,:]
        # subset and recorder adata_sel 
        adata_sel=adata[list(common_bc),:]
        # set name to pixel.x x pixel.y 
        pos_sel["idx"] = list(pos_sel[4].astype(str)  +"x" + pos_sel[5].astype(str))
        pos_sel.index = list(pos_sel["idx"])

        adata_sel.obs_names = list(pos_sel["idx"])
        mtx = pd.DataFrame(adata_sel.X.todense())
        mtx.index = adata_sel.obs_names
        mtx.columns = adata_sel.var_names
        
        print(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        mtx.to_csv("{}{}.csv".format(out_dir,f))
        pos_save = pos_sel.loc[:,[4,5]]
        pos_save.to_csv("{}{}_coor.csv".format(out_dir,f),header=False,index=False)
        
        
"""
    Helper function for loading organized 10X Visium output, loading result prepared by prepare_paste_input_dir()
    Input: 
        1. parent directory for all the inputs 
        2. sample names to be appended

    Assumption:
        1. Spot-gene count matrix is saved as data_dir + slice_name + ".csv"
        2. Coordinates are saved as data_dir + slice_name + "_coor.csv"

    Output: a list of AnnData objects, one for each slice 

    Example input: 
        data_dirs = ["../animals/animal14/"]*3 + ["../animals/animal53/"]*4
        slice_names = ["076D", "075B", "072B"] + ["106C","107C","108D","076A"]
        load_slices(data_dirs,slice_names)

"""

def load_slices(data_dirs, slice_names):
    
    slices = []  
    for i in range(len(slice_names)):
        data_dir = data_dirs[i]
        slice_name = slice_names[i]
        slice_i = sc.read_csv(data_dir + slice_name + ".csv")
        slice_i_coor = np.genfromtxt(data_dir + slice_name + "_coor.csv", delimiter = ',')
        slice_i.obsm['spatial'] = slice_i_coor
        # Preprocess slices
        sc.pp.filter_genes(slice_i, min_counts = 15)
        sc.pp.filter_cells(slice_i, min_counts = 100)
        slices.append(slice_i)
    return slices


def plot_expression(
    sliceX: AnnData, 
    gene: str, 
    title: str,
    ax: Optional[plt.Axes] = None,
    hue_norm = None,
    legend = False) -> None:
    """
    Plots slice spatial coordinates, colored by gene expression 
    
    Args:
        sliceX: Slice to be plotted.
        gene: gene name to be plotted
        title: title of the plot
    """
    exp = sliceX[:,gene].X.reshape(-1)
    color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
    pd_plot=pd.DataFrame(data={"exp":exp,
                           "x_cord":sliceX.obsm['spatial'][:,0],
                           "y_cord":sliceX.obsm['spatial'][:,1]})

    plot=sns.scatterplot(data=pd_plot,
                         x='y_cord',
                         y='x_cord',
                         hue='exp',
                         hue_norm = hue_norm,
                         #title=title,
                         palette=color_self,
                         size=5,
                         #size=100000/pd_plot.shape[0],
                         #legend="full",
                         linewidth=0,
                         ax=ax)
    plot.set(title=title)
    if not legend:
        plot.get_legend().remove()

    plot.axes.set_xlim([-1000, 1000])
    plot.axes.set_ylim([-1000, 1000])
    if ax:
        ax.invert_yaxis()
        #ax.axis('off')


def plot_multiple_expression(slices_to_plot,
                             plots_w,
                             plots_h, 
                             gene,
                             slice_names,
                             out_dir,
                             suffix="_combined",
                             save = True,
                             figsize=(22,22)):
    """
    Plots multiple slices' spatial coordinates , colored by gene expression 
    
    """
    
    plt.clf()
    w = plots_w
    h = plots_h
    exp_min = 1000
    exp_max = -1000

    # create subplots
    fig, axs = plt.subplots(w, h, figsize=figsize)
    # look for expression min and max among all plots to be plotted
    for i in range(len(slices_to_plot)):
        #print("(x,y) = ({},{})".format(str(x),str(y)))
        exp_i = slices_to_plot[i][:,gene].X.reshape(-1)
        exp_min = min(exp_min,np.min(exp_i))
        exp_max = max(exp_max,np.max(exp_i))
        #print("(min,max) = ({},{})".format(str(exp_min),str(exp_max)))

    # normalize the color based on min and max
    norm = plt.Normalize(exp_min, exp_max)
    color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
    sm = plt.cm.ScalarMappable(cmap=color_self, norm=norm)
    sm.set_array([])

    print("(min,max) = ({},{})".format(str(exp_min),str(exp_max)))

    # plot each slide as one suplot
    i = 0
    for x in range(w):
        for y in range(h):
            print("(x,y) = ({},{})".format(str(x),str(y)))
            if w == 1:
                ax_i = axs[y]
                # set global color bar position
                ax = axs[h-1]
            elif h == 1:
                ax_i = axs[x]
                ax = axs[w-1]
            else:
                ax_i = axs[x,y]
                ax = axs[0,h-1]
            if i < len(slices_to_plot):
                plot_expression(slices_to_plot[i],gene,slice_names[i],ax=ax_i,hue_norm=norm)
            else:
                fig.delaxes(ax_i)
            i +=1 

    fig.suptitle(gene, fontsize=16)

    # add color bar on the right, does not shrink the rightmost plot on the first row
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sm, cax=cax)

    if save: 
        # save 
        plt.savefig("{out}/expression_plot_{gene}{suffix}.png".format(out=out_dir,gene = gene,suffix=suffix, dpi=1200))
        plt.close()
    else:
        #plot 
        plt.show()












