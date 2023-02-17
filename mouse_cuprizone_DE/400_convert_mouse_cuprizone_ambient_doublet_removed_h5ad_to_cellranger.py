#conda activate /home/kli3/anaconda3/envs/cellxgeneVIP1.0.5/; python
import scanpy as sc
import pandas as pd
import os
import scipy.io

destination="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed"

#THIS DATASET HAS CELL TYPE LABELS PROCESSED BY SARBOTTAM PIYA ON 02/13/2023 
#UNKNOWN IF ANY QC WAS APPLIED
filename3="/edgehpc/dept/compbio/projects/TST11621/sarbottam/result/clustering/cuprizonemouse_annotated_Ambient_RNA_Doublet_removed.h5ad"
ad3=sc.read(filename3)
#The object has a raw.X and a metadata slot

#Write out full var slot for additional mapping columns
pd.DataFrame(ad3.var).to_csv(os.path.join(destination, "genes_full_info.tsv" ),   sep = "\t", index = True)
pd.DataFrame(ad3.obs.index).to_csv(os.path.join(destination, "barcodes.tsv"), sep = "\t", index = False)
ad3.obs.to_csv(os.path.join(destination, "metadata.tsv"), sep = "\t", index = True)
scipy.io.mmwrite(os.path.join(destination, "matrix.mtx"), ad3.raw.X) #access raw X for count

