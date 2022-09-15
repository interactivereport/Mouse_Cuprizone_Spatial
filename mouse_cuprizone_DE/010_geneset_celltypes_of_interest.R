#Load packages
library(tidyverse)
library(data.table)
library(glue)

set.seed(2022)

#User must set a program path
prg.dir <- file.path(getwd()) #"/edgehpc/dept/compbio/projects/SingleCell/programs/dev/mouse_cuprizone_DE_interactivereport/000_nebula_DE_func"

#NEBULA DE path
nebula_path <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE")
contrasts <- c("CRZ_vs_CTRL_LABELS","RCV_vs_CTRL_LABELS","CRZ_vs_RCV_LABELS")
ctois <- c("MFOL","MOL","NFOL","OPC")
ctois_regex <- paste0("(",paste(glue("{ctois}.csv"),collapse="|"),")")
nebula_de_files <- list.files(path=file.path(nebula_path, contrasts), pattern=ctois_regex, full.names = TRUE)
nebula_de_names <- gsub(".csv","",basename(nebula_de_files))

#Make an output directory for these
out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE/010_geneset_celltypes_of_interest")

source(file.path(prg.dir, "gene_set_func.R"))

for(i in 1:length(nebula_de_files)){
  
  neb_df <- fread(nebula_de_files[[i]], data.table=F)
  
  func_output_list <- gene_set_func(neb_df = neb_df,ora_fdr = 0.05,minGSSize = 10)
  
  fwrite(func_output_list[["ora"]], file.path(out_dir, glue("{nebula_de_names[[i]]}_ora.csv")))
  fwrite(func_output_list[["gsea"]], file.path(out_dir, glue("{nebula_de_names[[i]]}_gsea.csv")))
  
}
