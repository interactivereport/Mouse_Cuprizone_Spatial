import scanpy as sc
import pandas as pd
import os
import scipy.io

destination="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_LABELS"

#THIS DATASET HAS CELL TYPE LABELS PROCESSED BY MARK SHEEHAN ON 07/20/2022 
#THESE ARE THE FINAL LABELS
#THE CELLS WERE FILTERED
#FILTERS WERE >1000 GENES AND <10% MITOCHONDRIAL READS
filename3="/camhpc/home/msheehan/sc_sn_RNAseq/20220713_TST11621_cuprizone_finalLT/20220720_cuprizonemousebrain_cell2loc_and_OL_labels_added_clean.h5ad"
ad3=sc.read(filename3)

pd.DataFrame(ad3.var.index).to_csv(os.path.join(destination, "genes.tsv" ),   sep = "\t", index = False)
pd.DataFrame(ad3.obs.index).to_csv(os.path.join(destination, "barcodes.tsv"), sep = "\t", index = False)
ad3.obs.to_csv(os.path.join(destination, "metadata.tsv"), sep = "\t", index = True)
scipy.io.mmwrite(os.path.join(destination, "matrix.mtx"), ad3.raw.X) #access raw X for count

