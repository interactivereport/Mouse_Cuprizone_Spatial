#Load packages
library(tidyverse)
library(data.table)
library(glue)

#Set seed
set.seed(2023)

#Set directory
in_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed")

#Load metadata
meta_data <- fread(file.path(in_dir,"meta_format.tsv"))

#Tested cell types
cell_types_tested <- c("en_l6b","en_l5","en_hpc","in_sst","en","in_pvalb","en_thal","en_l56","in_vip","in_hypo","en_l6",
                       "in_lamp5","en_l23","en_l25","astrocytes","opc","nfol","endothelial","cop","microglia","mol","mfol")

#NEBULA DE path
nebula_path <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE")
contrasts <- paste0(c("CRZ_vs_CTRL","RCV_vs_CTRL","CRZ_vs_RCV"),"_ambient_doublet_removed")
ctois_regex <- paste0("(",paste(glue("{cell_types_tested}.csv"),collapse="|"),")")
nebula_de_files <- list.files(path=file.path(nebula_path, contrasts), pattern=ctois_regex, full.names = TRUE)
nebula_de_names <- gsub(".csv","",basename(nebula_de_files))

#Make an output directory for these
out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE/410_geneset_ambient_doublet_removed")

source("/edgehpc/dept/compbio/projects/SingleCell/programs/dev/mouse_cuprizone_DE/000_nebula_DE_func/gene_set_func.R")

library(parallel)

mclapply(nebula_de_files, function(x){
  neb_df <- fread(x, data.table=F)
  nebula_de_names <- gsub(".csv","",basename(x))
  func_output_list <- gene_set_func(neb_df = neb_df,ora_fdr = 0.05,minGSSize = 10)
  if(!is.null(func_output_list[["ora"]])){
    fwrite(func_output_list[["ora"]], file.path(out_dir, glue("{nebula_de_names}_ora.csv")))
  }
  if(!is.null(func_output_list[["gsea"]])){
    fwrite(func_output_list[["gsea"]], file.path(out_dir, glue("{nebula_de_names}_gsea.csv")))
  }
}, mc.cores = 4)
