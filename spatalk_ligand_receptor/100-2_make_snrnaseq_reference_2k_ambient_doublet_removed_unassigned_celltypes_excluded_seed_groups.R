library(tidyverse)
library(data.table)
library(spacexr)
library(Matrix)
library(doParallel)

# set random seed to make downsampled reference
# Create ref using Control, Cuprizone, and  Recovery snRNA-seq data
grp_arg=commandArgs(trailingOnly = T)[1]

# Single-Cell Reference
# In order to run RCTD, the first step is to process the single cell reference. 
# The reference is created using the RCTD Reference constructor function, which requires 3 parameters:
#   1. counts: A matrix (or dgCmatrix) representing Digital Gene Expression (DGE). 
#   Rownames should be genes and colnames represent barcodes/cell names. Counts should be untransformed count-level data.
# 2. cell_types: A named (by cell barcode) factor of cell type for each cell. The ‘levels’ of the factor would be the possible cell type identities.
# 3. nUMI: Optional, a named (by cell barcode) list of total counts or UMI’s appearing at each pixel. 
# If not provided, nUMI will be assumed to be the total counts appearing on each pixel.

counts <- readMM("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed/matrix.mtx")

meta_info=read.csv(file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed/meta_format_clean.csv"), stringsAsFactors=FALSE)
gene_info=fread(file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed/genes.tsv"), stringsAsFactors=FALSE, header=F) %>% #saved out as a headerless tsv w/2 columns
  data.frame()

colnames(counts) <- meta_info$cell
rownames(counts) <- gene_info[,1] 

# Filter meta data to one group
meta_info <- meta_info %>% 
  dplyr::filter(Group == grp_arg)

# Subset counts matrix
counts <- counts[, meta_info$cell]

#Reference() Takes dgCMatrix format or dense
counts_dgc <- as(counts, "dgCMatrix")
nUMI_obj <- setNames(nm = c("nUMI", "barcode"),stack(colSums(counts_dgc)))

rm(counts)

#Join the colsums of the spot object the metadata 
#Remove unassigned cell types from reference
meta_data <- meta_info %>% dplyr::select(barcode = cell, cluster = Celltype_annotation_clean) %>% left_join(nUMI_obj, by="barcode") %>%
  filter(!cluster %in% c("x", "not_assigned"))

#Make a named vector of cell types
cell_types <- meta_data$cluster
names(cell_types) <- meta_data$barcode # create cell_types named list

cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI
names(nUMI) <- meta_data$barcode # create nUMI named list

#Subset UMI data
counts_dgc <- counts_dgc[,meta_data$barcode, drop=FALSE]

### Create the Reference object
# downsample to 2000
set.seed(123)
reference <- Reference(counts = counts_dgc, 
                       cell_types = cell_types, 
                       nUMI = nUMI, 
                       n_max_cells = 2000, #Downsample to 2k max
                       min_UMI = 100 #Leave the default filter (should already be filtered)
)

## Save RDS object (optional)
saveRDS(reference, 
        file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/",
                  paste0('SCRef_2k_ambient_doublet_removed_unassigned_celltypes_excluded_seed_', grp_arg, '.rds')))
