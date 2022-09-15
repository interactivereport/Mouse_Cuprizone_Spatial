#Load packages
library(tidyverse)
library(data.table)
library(glue)

set.seed(2022)

#User must set a program path
prg.dir <- file.path(getwd()) #"/edgehpc/dept/compbio/projects/SingleCell/programs/dev/mouse_cuprizone_DE_interactivereport/000_nebula_DE_func"

metadata <- fread(file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/metadata_ST_region.tsv"))
annotation <- unique(metadata$annotation)

#NEBULA DE path
nebula_path <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE")
contrasts <- c("CRZ_vs_CTRL_ST","RCV_vs_CTRL_ST","CRZ_vs_RCV_ST")
annotation_use <- annotation[which(!annotation %in% "unassigned")]
ctois_regex <- paste0("(",paste(glue("{annotation_use}.csv"),collapse="|"),")")
nebula_de_files <- list.files(path=file.path(nebula_path, contrasts), pattern=ctois_regex, full.names = TRUE)
nebula_de_names <- gsub(".csv","",basename(nebula_de_files))

#Make an output directory for these
out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE/109_geneset_ST_annotation")

source(file.path(prg.dir, "gene_set_func.R"))

library(parallel)

mclapply(nebula_de_files, function(x){
  neb_df <- fread(x, data.table=F)
  nebula_de_names <- gsub(".csv","",basename(x))
  func_output_list <- gene_set_func(neb_df = neb_df,ora_fdr = 0.05,minGSSize = 10)
  fwrite(func_output_list[["ora"]], file.path(out_dir, glue("{nebula_de_names}_ora.csv")))
  fwrite(func_output_list[["gsea"]], file.path(out_dir, glue("{nebula_de_names}_gsea.csv")))
}, mc.cores = 4)
