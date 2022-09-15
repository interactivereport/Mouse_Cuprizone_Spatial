#Load packages
library(tidyverse)
library(data.table)
library(glue)

set.seed(2022)

#User must set a program path
prg.dir <- file.path(getwd()) #"/edgehpc/dept/compbio/projects/SingleCell/programs/dev/mouse_cuprizone_DE_interactivereport/000_nebula_DE_func"

metadata <- fread(file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_LABELS/metadata.tsv"))
final_label <- unique(metadata$final_label)

final_label_tested <- c("Astro","COP","Endo","Ext_ClauPyr",
                        "Ext_Hpc_CA1","Ext_Hpc_CA2","Ext_Hpc_CA3","Ext_Hpc_DG1","Ext_Hpc_DG2","Ext_L23",
                        "Ext_L25","Ext_L5_1","Ext_L5_2","Ext_L56","Ext_L6","Ext_L6B","Ext_Med","Ext_Pir",
                        "Ext_Thal_1","Ext_Thal_2","Inh_2","Inh_3","Inh_4","Inh_Lamp5","Inh_Meis2_2","Inh_Pvalb","Inh_Sst","Inh_Vip",
                        "MFOL","Micro","MOL","NFOL","OPC")

#NEBULA DE path
nebula_path <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE")
contrasts <- c("CRZ_vs_CTRL_LABELS","RCV_vs_CTRL_LABELS","CRZ_vs_RCV_LABELS")
ctois <- c("MFOL","MOL","NFOL","OPC")
final_label_use <- final_label_tested[which(!final_label_tested %in% ctois)]
ctois_regex <- paste0("(",paste(glue("{final_label_use}.csv"),collapse="|"),")")
nebula_de_files <- list.files(path=file.path(nebula_path, contrasts), pattern=ctois_regex, full.names = TRUE)
nebula_de_names <- gsub(".csv","",basename(nebula_de_files))

#Make an output directory for these
out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE/012_geneset_final_label")

source(file.path(prg.dir, "gene_set_func.R"))

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

#Get labels not tested
final_label_untested <- setdiff(final_label, final_label_tested)

#Do the test count to justify the removal of the cell types that were not tested
min_celltype_count <- metadata %>% 
  mutate(num_id = paste0(Number,"_",Group)) %>%
  count(num_id,final_label) %>% 
  complete(num_id, final_label, fill=list(n=0)) %>% #dplyr count and group_by drop empty factor levels
  filter(n < 10) %>% 
  pull(final_label) %>% unique()
min_celltype_count <- c("Unknown", min_celltype_count) #Unknown was dropped by default

all.equal(sort(final_label_untested), sort(min_celltype_count))

message(glue("Did not perform DE on {paste(min_celltype_count,collapse=',')}."))
