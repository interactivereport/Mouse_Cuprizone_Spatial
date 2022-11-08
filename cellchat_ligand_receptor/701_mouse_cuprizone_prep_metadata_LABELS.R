require(tidyverse)
require(magrittr)
library(data.table)
library(Matrix)

#Set the directory
data_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_LABELS")
out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/400_lr_pipeline")
#Load the metadata, barcodes, gene info for mouse cuprizone
meta_data <- fread(file.path(data_dir, "metadata.tsv"))

#Cell column out of index
#Remove special characters from metadata columns
meta_data_format <- meta_data %>%
  dplyr::rename(cell = V1) 

#Write out
fwrite(meta_data_format, file.path(out_dir, "701_mouse_cuprizone_meta_format.tsv"))

#Format feature file
gene_info <- fread(file.path(data_dir, "genes.tsv"))
gene_info$dummy_ensembl_id <- gene_info$index
write.table(gene_info, file.path(data_dir,"genes.tsv"), col.names=F, row.names=F, sep="\t", quote=F)

#Format barcodes
barcodes <- fread(file.path(data_dir, "barcodes.tsv"))
barcodes <- barcodes %>% filter(V1 != "0") #Formatting issue during extraction from h5ad
write.table(barcodes, file.path(data_dir,"barcodes.tsv"), col.names=F, row.names=F, sep="\t", quote=F)

#Format matrix
counts <- readMM(file.path(data_dir, "matrix.mtx"))
counts <- t(counts)
writeMM(counts, file.path(data_dir, "matrix.mtx"))