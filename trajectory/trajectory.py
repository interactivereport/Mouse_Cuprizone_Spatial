import scanpy as sc
import numpy as np
import gc
import pandas as pd

import scFates as scf

adata = sc.read_h5ad('/edgehpc/dept/compbio/projects/TST11621/sarbottam/result/Oligodendrocytes_subclustered.h5ad')
adata.X = adata.raw.X.copy()
adata = adata[adata.obs.Celltype_annotation=='MFOL']
sc.pp.filter_genes(adata,min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000,  min_mean=0.0125, max_mean=3, min_disp=0.5)
adata=adata[:,adata.var.highly_variable]
n_comps = 50
sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)


scf.tl.curve(adata,Nodes=50, use_rep="X_pca", ndims_rep=n_comps)
adata.obs['MFOL_Group'] = adata.obs['Group']
scf.pl.graph(adata, basis="pca", color_cells='MFOL_Group'), save='_trajectory_MFOL_DiseaseGroups.png')
