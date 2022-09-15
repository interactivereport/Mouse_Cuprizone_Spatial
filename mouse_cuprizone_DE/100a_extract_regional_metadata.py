import scanpy as sc
import pandas as pd
import os
import scipy.io

destination="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st"

#The anndata file contains metadata in obs with annotation information for each spot
filename3="/camhpc/ngs/projects/TST11523/notebook/TST11523_cVIP_3_customized_HHT-TUMAOYKH-general.h5ad"
ad3=sc.read(filename3)

ad3.obs.to_csv(os.path.join(destination, "metadata_ST_region.tsv"), sep = "\t", index = True)
