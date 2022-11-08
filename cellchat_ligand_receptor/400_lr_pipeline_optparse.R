#!/usr/bin/env Rscript

#Clear
rm(list=ls())

#Libraries
suppressPackageStartupMessages({
library(tidyverse)
library(magrittr)
library(Seurat)
library(data.table)
library(ggrepel)
library(scDiffCom)
library(optparse)
library(CellChat)
library(patchwork)
})

#Option list
option_list = list(

  #Set up input
  make_option("--input_10X", action = "store", default=NA, type="character",
  help="path to cellranger v2 input.  Requires matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files provided by 10X."),
  make_option("--input_meta", action = "store", default=NA, type="character",
  help="path to metadata.  Must have column 'cell' that matches the UMI count matrix."),
  
  #Set up contrast
  make_option("--contrast_col", action="store", default=NA, type="character",
  help="name of contrast column.  Must exist in metadata."),
  make_option("--contrast_ref", action="store", default=NA, type="character",
  help="name of contrast reference group.  Must exist in --contrast_col."),
  make_option("--contrast_alt", action="store", default=NA, type="character",
  help="name of contrast alternative group.  Must exist in --contrast_col."),
  make_option("--group_col", action="store", default=NA, type="character",
  help="name of grouping column.  Must exist in metadata."),
  
  #Tool selection and analysis settings
  make_option("--input_species", action = "store", default=NA, type="character",
  help="Species for input data.  For CellChat or scDiffCom, select mouse or human.  scDiffCom allows rat.  CellChat allows zebrafish."),
  
  make_option("--cellchat", action = "store_true", default=FALSE, type="logical",
  help="include this flag to run CellChat"),
  make_option("--scdiffcom", action = "store_true", default=FALSE, type="logical",
  help="include this flag to run scDiffCom"),
    make_option("--scdiffcom_parallel", action="store_true", default=FALSE, type="logical",
    help="include this flag for scDiffCom parallel processing for larger datasets."),
  
  #File output
  make_option("--out_dir", action = "store", default=NA, type="character",
  help="path to output directory for LR pipeline")
   
)

opt = parse_args(OptionParser(option_list=option_list))

############# Parse input #############
#read in metadata
meta_data=fread(file.path(opt$input_meta))

#Read in counts
dat <- Read10X(file.path(opt$input_10X))

############# scDiffCom ###############

if(opt$scdiffcom){

  #Format the metadata
  meta_format <- meta_data %>% 
    dplyr::filter(!!sym(opt$contrast_col) %in% c(opt$contrast_ref, opt$contrast_alt)) 
  
  #Create the matrix
  order_index = which(colnames(dat) %in% meta_format$cell)
  order_counts = dat[, order_index]
  
  #Create the matrix
  order_index = which(colnames(dat) %in% meta_format$cell)
  order_counts = dat[, order_index]
  rm(dat)
    
  #Make the Seurat
  meta_info_format <- meta_format %>% column_to_rownames("cell")
  srt_obj <- CreateSeuratObject(counts = order_counts, meta.data = meta_info_format)
  
  #Process data (default handling will back-convert from norm. log1p in data slot of RNA assay)
  srt_obj <- NormalizeData(srt_obj, assay="RNA") 
  
  if(opt$scdiffcom_parallel){
    library(future)
    plan(multicore, workers=4) #match nodes
    options(future.globals.maxSize= 891289600) #Set for 850mb 
  }
  
  #Do default handling of the data
  scdiffcom_obj <- run_interaction_analysis(
    seurat_object = srt_obj,
    LRI_species = opt$input_species, #Either "mouse", "human" or "rat"
    seurat_celltype_id = opt$group_col, 
    seurat_condition_id = list(
      column_name = opt$contrast_col,
      cond1_name = opt$contrast_ref,
      cond2_name = opt$contrast_alt
    )
  )
  
  #Save the scdiffcom object
  saveRDS(scdiffcom_obj, file.path(opt$out_dir, "scdiffcom_obj.rds"))
  
  #Save the detected CCI table
  CCI_detected <- GetTableCCI(scdiffcom_obj, type = "detected", simplified = TRUE)
  fwrite(CCI_detected, file.path(opt$out_dir, "scdiffcom_CCI_detected.csv"))
  
  #Do the ORA
  ORA_results <- GetTableORA(scdiffcom_obj, categories = "all", simplified = TRUE)
  
  #Save the LR ORA table
  fwrite(ORA_results$LRI, file.path(opt$out_dir, "scdiffcom_ORA_LRI.csv"))
  #Save the ER ORA table
  fwrite(ORA_results$ER_CELLTYPES, file.path(opt$out_dir, "scdiffcom_ORA_ER.csv")) 
  #Save the GO ORA table
  fwrite(ORA_results$GO_TERMS, file.path(opt$out_dir, "scdiffcom_ORA_GO.csv"))  

}

############# CellChat ###############

if(opt$cellchat){

  #CellChat does per-group analysis 
  meta_ref <- meta_data %>%
    dplyr::filter(!!sym(opt$contrast_col) %in% c(opt$contrast_ref)) #Select reference group of interest
    
  meta_alt <- meta_data %>%
    dplyr::filter(!!sym(opt$contrast_col) %in% c(opt$contrast_alt)) #Select alternative group of interest

  #Format the counts
  index_ref = which(colnames(dat) %in% meta_ref$cell)
  filter_counts_ref = dat[, index_ref]

  index_alt <- which(colnames(dat) %in% meta_alt$cell)
  filter_counts_alt = dat[, index_alt]

  rm(dat)

  #Make the Seurat objects
  meta_info_ref_format <- meta_ref %>% column_to_rownames("cell")
  meta_info_alt_format <- meta_alt %>% column_to_rownames("cell")
  srt_ref <- CreateSeuratObject(counts = filter_counts_ref, meta.data = meta_info_ref_format)
  srt_alt <- CreateSeuratObject(counts = filter_counts_alt, meta.data = meta_info_alt_format)

  #Process data 
  srt_ref <- NormalizeData(srt_ref, assay="RNA")
  srt_alt <- NormalizeData(srt_alt, assay="RNA")
  srt_ref <- FindVariableFeatures(srt_ref, nfeatures=2000)
  srt_alt <- FindVariableFeatures(srt_alt, nfeatures=2000)

  #Set ident
  Idents(srt_ref) <- meta_ref[, opt$group_col]
  Idents(srt_alt) <- meta_alt[, opt$group_col]
  
  #Create a function to do the rest of the paired steps
  cellchat_fun <- function(srt_obj, species, min_num_cells){
    #Create the object
    cellchat_obj <- createCellChat(object = srt_obj,
                                     group.by = "ident",
                                     assay = "RNA")                                  
    #Set the database
    if(species == "human"){                                  
      CellChatDB <- CellChatDB.human 
    }
    if(species == "mouse"){
      CellChatDB <- CellChatDB.mouse 
    }
    if(species == "zebrafish"){
      CellChatDB <- CellChatDB.zebrafish
    }
    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling, hardcoded 
    
    # set the used database in the object
    cellchat_obj@DB <- CellChatDB.use
    
    #Do the analytical steps
    cellchat_obj <- subsetData(cellchat_obj)
    cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
    cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
    # change default type of triMean to truncatedMean to have more interactions on Aug 24, 2022
    # set population.size=True on Sep 26, 2022
    cellchat_obj <- computeCommunProb(cellchat_obj, type = "truncatedMean", trim = 0.1, population.size = TRUE) 
    cellchat_obj <- filterCommunication(cellchat_obj, min.cells = min_num_cells)
    cellchat_obj <- computeCommunProbPathway(cellchat_obj)
    cellchat_obj <- aggregateNet(cellchat_obj)
  
  }
  
  #Run the wrapper function on each contrast
  cellchat_ref <- cellchat_fun(srt_ref, species=opt$input_species, min_num_cells = 10)
  cellchat_alt <- cellchat_fun(srt_alt, species=opt$input_species, min_num_cells = 10)
  
  #Merge the processed contrast objects
  object.list <- list(ref = cellchat_ref, alt = cellchat_alt)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))

  # deal with the situation where datasets have different cell type compositions
  if (length(levels(cellchat_alt@idents)) != length(levels(cellchat_ref@idents))) {
    cellchat <- liftCellChat(cellchat)
  }
  
  #Write out the merged cell chat object
  saveRDS(cellchat_ref, file.path(opt$out_dir, "cellchat_ref.rds"))
  saveRDS(cellchat_alt, file.path(opt$out_dir, "cellchat_alt.rds"))
  saveRDS(cellchat, file.path(opt$out_dir, "cellchat.rds"))  
  
  #Network comparison table (extract net visual)
  netvisual_use <- unique(cellchat@idents$joint)
  
  gg1 <- netVisual_bubble(cellchat, sources.use = netvisual_use,
                        targets.use = netvisual_use,  
                        comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased in alternative group",
                        angle.x = 45, remove.isolate = T,
                        return.data=T)

  gg2 <- netVisual_bubble(cellchat, sources.use = netvisual_use,
                          targets.use = netvisual_use,  
                          comparison = c(1, 2), max.dataset = 1, 
                          title.name = "Decreased in alternative group", 
                          angle.x = 45, remove.isolate = T,
                          return.data=T)

  
  # save dataframe
  df1 <- gg1[["communication"]] %>% 
  data.frame() %>% 
  mutate(Regulated = "up")

  df2 <- gg2[["communication"]] %>% 
    data.frame() %>% 
    mutate(Regulated = "down") %>% 
    rbind(df1)

  fwrite(df2, file.path(opt$out_dir, "cellchat_compare_interaction_dif_signaling.csv"))
  
}
