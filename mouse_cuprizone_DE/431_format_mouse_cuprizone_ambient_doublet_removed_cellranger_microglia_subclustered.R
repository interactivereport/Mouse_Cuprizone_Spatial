require(tidyverse)
require(magrittr)
library(data.table)
library(Matrix)

#Set the directory
data_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed/microglia_subclustered")
out_dir <- data_dir

#Load the metadata, barcodes, gene info for mouse cuprizone
meta_data <- fread(file.path(data_dir, "metadata.tsv"))

#Cell column out of index
#Remove special characters from metadata columns
#Format additional clustering groups
meta_data_format <- meta_data %>%
  dplyr::rename(cell = V1) %>%
  dplyr::mutate(leiden_contrast = paste0("cluster_",leiden),
                leiden_contrast3 = case_when(leiden %in% c(1,2) ~ "cluster_1_2",
                                             leiden %in% c(0) ~ "cluster_0",
                                             TRUE ~ "not_used"),
                leiden_contrast4 = case_when(leiden %in% c(5) ~ "cluster_5",
                                             TRUE ~ "rest")
                )

#Write out
fwrite(meta_data_format, file.path(out_dir, "meta_format.tsv"))

#Format feature file
gene_info <- fread(file.path(data_dir, "genes_full_info.tsv"))
gene_info <- gene_info %>% dplyr::select(index = V1)
gene_info$dummy_ensembl_id <- gene_info$index
write.table(gene_info, file.path(data_dir,"genes.tsv"), col.names=F, row.names=F, sep="\t", quote=F)

#Format barcodes
barcodes <- fread(file.path(data_dir, "barcodes.tsv"))
barcodes <- barcodes %>% filter(V1 != "0") 
write.table(barcodes, file.path(data_dir,"barcodes.tsv"), col.names=F, row.names=F, sep="\t", quote=F)

#Format matrix
counts <- readMM(file.path(data_dir, "matrix.mtx"))
counts <- t(counts)
writeMM(counts, file.path(data_dir, "matrix.mtx"))
