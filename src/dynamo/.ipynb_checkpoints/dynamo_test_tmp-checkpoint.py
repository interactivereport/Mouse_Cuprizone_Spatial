import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import dynamo as dyn
from dynamo.dynamo_logger import main_info, LoggerManager
import gc
import warnings

warnings.filterwarnings('ignore')
LoggerManager.main_logger.setLevel(LoggerManager.INFO)
dyn.configuration.set_figure_params('dynamo', background='white') # jupter notebooks
# dyn.configuration.set_figure_params('dynamo', background='black') # presentation
dyn.get_all_dependencies_version()

# def hhhh():
#     filename = '/edgehpc/dept/compbio/users/msheehan/msheehan/sc_sn_RNAseq/20210601_TST11621_cuprizone_mouse_brain_snRNAseq/analysis.01.clustering/res.04.clustering_sub1234/cuprizonemouse_brain_alllevels_comb_v3_subclusterglia_clean.h5ad'

def hh():
    adata = dyn.read('/edgehpc/dept/compbio/users/whu1/spatial/data/cpz_intro_exon.h5ad')
    adata.obs['Cell_type'] = adata.obs['final_label'] 
    adata.layers['spliced'] = adata.layers['intronic']
    adata.layers['unspliced'] = adata.layers['exonic']
    del adata.layers["intronic"]
    del adata.layers["exonic"]
    gc.collect()
    print(adata)
    dyn.pp.recipe_monocle(adata, keep_filtered_cells=False, keep_filtered_genes=True, maintain_n_top_genes=True)
    print(adata)
    
    dyn.tl.dynamics(adata, model='deterministic', cores=400)
    # dyn.pl.umap(adata, color='Cell_type')
    
    # dyn.pl.cell_wise_vectors(adata, color=['Cell_type'], basis='umap', show_legend='on data', quiver_length=6, quiver_size=6, pointsize=0.1, show_arrowed_spines=False, save_show_or_return='save', s_kwargs_dict={'dpi':300, 'figsize':(12,8)})
    
    dyn.pl.streamline_plot(adata, color=['Cell_type'], basis='umap', show_legend='on data', show_arrowed_spines=True, save_show_or_return='save', s_kwargs_dict={'dpi':300, 'figsize':(12,8)})
    

if __name__ == '__main__':
    hh()
