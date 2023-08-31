library(tidyverse)
library(data.table)
library(Matrix)
library(doParallel)
library(ggplot2)
library(Seurat)
library(glue)
library(ggpubr)
library(ggrepel)
library(BayesSpace, lib.loc="/home/ezhao1/R/x86_64-pc-linux-gnu-library/4.1")
library(spacexr, lib.loc = "/home/jluo/R/x86_64-pc-linux-gnu-library/4.2") # updated version 2.2.0
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  # Select SVG or DEG analysis
  make_option("--SVG", action = "store_true", default = FALSE, type = "logical",
              help = "Include this flag to perform SVG analysis [required]"),
  make_option("--DEG", action = "store_true", default = FALSE, type = "logical",
              help = "Include this flag to perform DEG analysis [required]"),
  
  # Select regions, groups, and cell types
  make_option("--group", action = "store", default = NA, type = "character",
              help = "Specify the group in SVG analysis [required]"),
  make_option("--group_ref", action = "store", default = NA, type = "character",
              help = "Specify the ref group in DEG analysis [required]"),
  make_option("--group_alt", action = "store", default = NA, type = "character",
              help = "Specify the alt group in DEG analysis [required]"),
  make_option("--region_ref", action = "store", default = NA, type = "character",
              help = "Specify the ref region [required]"),
  make_option("--region_alt", action = "store", default = NA, type = "character",
              help = "Specify the alt region [required]"),
  make_option("--region", action = "store", default = NA, type = "character",
              help = "Specify the region in DEG analysis [required]"),
  make_option("--celltypes", action = "store", default = NA, type = "character",
              help = "Specify the cell types of interest [required]"),
  make_option("--samp_excl", action = "store", default = NA, type = "character",
              help = "Specify the samples to exclude because of poor quality [required]"),
  
  # Parameters
  make_option("--log_fc_thresh", action = "store", default = 0.4, type = "numeric",
              help = "Log fold change threshold for population inference [required]"),
  make_option("--cell_type_count_thresh", action = "store", default = 0, type = "numeric",
              help = "cell type count threshold [required]"),
  make_option("--weight_thresh", action = "store", default = 0, type = "numeric",
              help = "Minimal sum of cell type weight in spots [required]"),
  make_option("--gene_threshold", action = "store", default = 0.00005, type = "numeric",
              help = "Minimal expression of genes [required]"),
  make_option("--max_core", action = "store", default = 4, type = "integer",
              help = "Number of cores in parallel computing [required]"),
  
  # Input
  make_option("--in_dir", action = "store", default = NA, type = "character",
              help = "The directory of input files [required]"),
  make_option("--ct_wt", action = "store", default = NA, type = "character",
              help = "The deconvoluted cell type weight matrix [required]"),
  make_option("--sc_ref", action = "store", default = NA, type = "character",
              help = "The directory of single cell reference file [required]"),
  make_option("--puck_dir", action = "store", default = NA, type = "character",
              help = "The directory of puck files [required]"),
  make_option("--qc_spots", action = "store", default = NA, type = "character",
              help = "The QCed spots [required]"),
  make_option("--qc_genes", action = "store", default = NA, type = "character",
              help = "The QCed genes [required]"),
  make_option("--meta_data", action = "store", default = NA, type = "character",
              help = "The meta data [required]"),
  
  # Population inference approach selection
  make_option("--meta_regression", action = "store_true", default = FALSE, type = "logical",
              help = "Using meta regression approach for population inference [required]"),
  make_option("--ave_pop", action = "store_true", default = FALSE, type = "logical",
              help = "Using average approach for population inference [required]"),
  
  # Output
  make_option("--out_dir", action = "store", default = NA, type = "character",
              help = "The output directory [required]")
)))


# Turn off cell type filter and weight_thresh filter to not exclude spots and cell types
# Some parameters
log_fc_thresh = opt$log_fc_thresh
cell_type_count_thresh = opt$cell_type_count_thresh
weight_thresh = opt$weight_thresh
gene_threshold = opt$gene_threshold
max_core=opt$max_core

# Specify cell types to analyze
cell_types <- unlist(strsplit(opt$celltypes, split="[ ,;]+"))

# Debug
if(FALSE){
  group="CPZ_CTL"
  region=c("TH","hypoTH")
  opt=list(group="CTL",
           region_ref="CC", region_alt="cortex",
           # region="cortex",
           celltypes="Ext_Hpc_CA1,Ext_L5_1,Astro,OPC,Ext_L23,MFOL,Micro,MOL",
           sc_ref="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/SCRef_1k_ambient_doublet_removed.rds",
           puck_dir="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS/spacexr_format",
           qc_spots="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/443_ST_QC_summary/443_ST_QC_scenario4_spots.csv",
           qc_genes="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/443_ST_QC_summary/443_ST_QC_scenario4_genes.csv",
           meta_data="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS/metadata.tsv",
           ct_wt="/edgehpc/dept/compbio/projects/TST11523/2023/cell2loc/result/cell_type_weight_mtx_all.csv",
           samp_excl="072A,072D,074D",
           group_ref="CTL",
           group_alt="CPZ",
           region="TH")
}

################################
#### SVG analysis ##############
################################

if(opt$SVG){
  
  # Specify the regions of contrast in a group
  grps <- opt$group
  region=c(opt$region_alt, opt$region_ref)
  
  # The file prefix for outputs
  nm <- paste(opt$group, opt$region_alt, opt$region_ref, sep = "_")

  
  out_dir=file.path(opt$out_dir, nm)
  
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
    system(glue::glue("chmod 775 {out_dir}"))
  }
  
  st_meta_tmp <- fread(opt$meta_data) %>%
    dplyr::select(sample, type) %>% distinct() %>%
    dplyr::filter(type %in% grps)

  # Load the sc ref
  sc_ref <- readRDS(file.path(opt$sc_ref))
  
  # Directory of puck files created
  puck_dir <- opt$puck_dir
  
  # QCed barcodes and genes
  qc_spots <- fread(file.path(opt$qc_spots))
  qc_genes <- fread(file.path(opt$qc_genes))
  
  
  #Function that will prepare the correct files to use
  puck_prep <- function(puck_dir, sampleid_vector, annot, do_qc_filter=F, qc_barcodes=NULL, qc_genes=NULL){
    #Prepare lookup regex
    sampleid_vector_regex <- paste0("(",paste(sampleid_vector, collapse="|"),")",glue("_{annot}_puck.rds$"))
  
    #Source the spacexr puck format
    puck_file_list <- list.files(path = puck_dir, pattern=sampleid_vector_regex, full.names=T)
  
    #Read in found puck files
    puck_list <- lapply(puck_file_list, readRDS)
  
    names(puck_list) <- gsub("_puck.rds","",basename(puck_file_list))
  
    if(do_qc_filter){
      stopifnot( !is.null(qc_barcodes) & !is.null(qc_genes) )
      #Spot-level
      puck_list <- lapply(puck_list, function(x){
        spacexr::restrict_puck(puck = x,
                               barcodes = qc_barcodes)})
      
      #Gene-level
      puck_list <- lapply(puck_list, function(x){
        # modified Feb 9, 2023
        # make sure qc_genes are in the puck file, otherwise,can not subset the matrix with index not exists.
        qc_gene_present <- qc_genes[qc_genes %in% x@counts@Dimnames[[1]]]
        spacexr::restrict_counts(puck = x,
                                 gene_list = qc_gene_present, UMI_thresh = 0, UMI_max = Inf, counts_thresh = 0)})
    }
    return(puck_list)
  }
  # End function: puck_prep
  
  # Get and process puck files
  puck_list <- puck_prep(puck_dir = puck_dir,
                         sampleid_vector = st_meta_tmp %>% pull("sample") %>% unique(),
                         annot = paste(region, collapse = "_"),
                         qc_barcodes = qc_spots$sampleid_spot,
                         qc_genes = qc_genes$gene_name,
                         do_qc_filter=T)
  
  
  replicate_names <- str_extract(names(puck_list), ".*?(?=_)")
  
  group_ids <- as.numeric(factor(st_meta_tmp$type[which(st_meta_tmp$sample %in% replicate_names)], levels = c("CTL","RCV","CPZ"))) #group ids must be numeric
  
  # Create RCTD object  
  myRCTD.reps <- create.RCTD.replicates(spatialRNA.replicates = puck_list,
                                        reference = sc_ref,
                                        replicate_names = replicate_names,
                                        group_ids = group_ids,
                                        max_cores = max_core,
                                        keep_reference = T,
                                        CELL_MIN_INSTANCE = 0 #Don't drop celltypes on the basis of minimum number
                                        )
  
  saveRDS(myRCTD.reps, file.path(out_dir, paste0(nm,"_rctd_create_replicates_full_no_filter.rds")))
  
  # Using customer cell type weight matrix
  c2l_wts <- fread(file.path(opt$ct_wt)) %>%
    column_to_rownames("V1") 

  # Replace NAs in weights with 0
  c2l_wts[is.na(c2l_wts)] <- 0
  
  # Scale weight matrix to have row sum of 1
  c2l_wts_scaled_mtx <- sweep(c2l_wts, 1, rowSums(c2l_wts), '/')
  
  RCTD_controls <- myRCTD.reps
  
  # Import weight for each replicate
  for (i in 1:length(RCTD_controls@RCTD.reps)) {
    # Subset each replicate barcodes
    barcodes <- RCTD_controls@RCTD.reps[[i]]@spatialRNA@counts@Dimnames[[2]]
    c2l_wt_rep <- c2l_wts_scaled_mtx %>%
      dplyr::filter(rownames(.) %in% barcodes)
  
    RCTD_controls@RCTD.reps[[i]] <- import_weights(RCTD_controls@RCTD.reps[[i]], c2l_wt_rep)
  
    # Need to set the doublet_mode as "full" so that the code will use this weights in slot: myRCTD@results$weights
    RCTD_controls@RCTD.reps[[i]]@config$doublet_mode <- "full"
    RCTD_controls@RCTD.reps[[i]]@config$RCTDmode <- "full"
  
  }
  

  # Create exvar file
  puck_exvar <- list()

  st_meta_harmonized <- read.csv("/mnt/depts/dept04/compbio/projects/SingleCell/results/dev/spatial_tx/cside/st_meta_harmonized.csv") %>%
    dplyr::filter(annotation %in% region) %>%
    dplyr::filter(type %in% grps)
  
  
  for(i in 1:length(myRCTD.reps@RCTD.reps)){
    barcodes <- myRCTD.reps@RCTD.reps[[i]]@spatialRNA@counts@Dimnames[[2]]
    exvar <- as.integer(st_meta_harmonized[st_meta_harmonized$sampleid_spot %in% barcodes,]$annotation == region[1])
    names(exvar) <- barcodes
    puck_exvar[[i]] <-  exvar
  }
  
  ## run cside  
  myRCTD.reps <- run.CSIDE.replicates(RCTD_controls,
                                      cell_types =  cell_types,
                                      puck_exvar,
                                      weight_threshold = weight_thresh,
                                      log_fc_thresh = log_fc_thresh,
                                      doublet_mode = FALSE, #Ran with full weights, not doublet
                                      gene_threshold = gene_threshold,
                                      population_de = F, # Get DE across all cell types, run this later
                                      cell_type_threshold = cell_type_count_thresh,
                                      fdr = 0.05)
  
  saveRDS(myRCTD.reps,
          file.path(out_dir, paste0(nm,"_cside_full_no_filter_no_popde_c2l.rds")))
  
  # Use average method for population inference
  if(opt$ave_pop){
    myRCTD.reps <- CSIDE.population.inference(myRCTD.reps, 
                                              log_fc_thresh=log_fc_thresh)

    saveRDS(myRCTD.reps, file.path(out_dir, paste0(nm,"_cside_full_no_filter_popde_c2l_average.rds")))
    
    res_df <- map(myRCTD.reps@population_de_results, function(x) x %>% rownames_to_column("gene")) %>%
      rbindlist(idcol="cell") %>%
      mutate(log10p=-log10(p),
             log2FC= log2(exp(log_fc_est))) 
    
    write.csv(res_df, file.path(out_dir,
                                paste0(nm,"_SVG_full_no_filter_popde_c2l_average_pop.csv")))
  }
  
  # Use meta regression for population inference
  # Note: this method is not ready in the current version of Spacexr package
  if(opt$meta_regression){
    ## Run meta regression to compare groups
    
    meta_design_matrix <- st_meta_tmp %>%
      dplyr::filter(sample %in% names(myRCTD.reps@group_ids)) %>% 
      mutate(type=factor(type, levels=c(grps[[2]], grps[[1]]))) %>% 
      column_to_rownames("sample") %>% 
      model.matrix(~type,.) 
      
     
    myRCTD.reps <- CSIDE.population.inference(myRCTD.reps, log_fc_thresh=log_fc_thresh,
                                              meta = T,
                                              meta.design.matrix = meta_design_matrix,
                                              fdr = 0.05,
                                              meta.test_var = paste0("type", grps[[1]]))
    
    
    
    saveRDS(myRCTD.reps, file.path(out_dir, paste0(nm,"_cside_full_no_filter_popde_c2l_meta_regression.rds")))
    
    res_df <- map(myRCTD.reps@population_de_results, function(x) x %>% rownames_to_column("gene")) %>%
      rbindlist(idcol="cell") %>%
      mutate(log10p=-log10(p),
             log2FC= log2(exp(log_fc_est))) %>%
      mutate(note=case_when(log2FC > 0.5 & p > -log10(0.01) ~ "up",
                            log2FC < -0.5 & p > -log10(0.01) ~ "down",
                            TRUE ~ "not sig"),
             note=factor(note, levels = c("up","down","not sig")))
    
    write.csv(res_df, file.path(out_dir,
                                paste0(nm,"_SVG_full_no_filter_popde_c2l_meta_regression.csv")))    
    
  }
}


################################
#### DEG analysis ##############
################################
if(opt$DEG){
  
  ## run DEG analysis using pooled data with animal_id as covariate 
  # Select groups of contrast in a region
  grps <- c(opt$group_alt, opt$group_ref)
  region=opt$region
  
  # File prefix for outputs
  nm <- paste(opt$group_alt, opt$group_ref, opt$region, sep = "_")

  out_dir=file.path(opt$out_dir, nm)
  
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }
  
  
  # Load reference
  reference <- readRDS(file.path(opt$sc_ref))
  
  # Create puck files
  set.seed(100)
  ID.map = read.csv("/edgehpc/dept/compbio/users/ezhao/data-raw/ID.map.csv")
  ID.map$slide = substr(ID.map$SampleID, 1,3) #get slide number
  ID.map$slide.letter = substr(ID.map$SampleID, 4,4)
  method = "harmony"
  
  # Select samples in groups
  ID.map <- ID.map %>%
    dplyr::filter(Type %in% grps)
  
  # exclude more samples if --samp_excl is not NA
  if(! is.na(opt$samp_excl)){
    ID.map <- ID.map %>% 
      dplyr::filter(! SampleID %in% unlist(strsplit(opt$samp_excl, split=",| ")))
  }

  st_meta_harmonized <- read.csv("/mnt/depts/dept04/compbio/projects/SingleCell/results/dev/spatial_tx/cside/st_meta_harmonized.csv")
  st_meta_harmonized <- st_meta_harmonized %>%
    dplyr::filter(annotation == region, type %in% grps)
  
  # exclude samples if --samp_excl is not NA
  if(! is.na(opt$samp_excl)){
    st_meta_harmonized <- st_meta_harmonized %>%
      mutate(SampleID=str_extract(sampleid_spot, ".*?(?=_)")) %>% 
      dplyr::filter(! SampleID %in% unlist(strsplit(opt$samp_excl, split=",| "))) 
  }
  
  spot_use <- st_meta_harmonized %>%
    pull(sampleid_spot)


  # create one puck file
  sample_meta_all <- data.frame()
  n_gene <- readVisium(paste0("/edgehpc/dept/compbio/users/ezhao/data-raw/", ID.map$SampleID[1]))@assays@data$counts@Dim[1]
  
  sample_count_all <- data.frame(matrix(nrow = n_gene))
  
  for (i in 1:nrow(ID.map)){
    sample.id = ID.map$SampleID[i]
    sample.sce = readVisium(paste0("/edgehpc/dept/compbio/users/ezhao/data-raw/", sample.id))
    
    # Extract raw counts and metadata, format the index
    sample_metadata <- data.frame(colData(sample.sce),stringsAsFactors = FALSE) %>%
      mutate( sampleid_spot = paste0(sample.id,"_",spot) ) %>%
      dplyr::filter(sampleid_spot %in% spot_use)
    
    sample_meta_all <- rbind(sample_meta_all, sample_metadata)
    
    sample_raw_count <- sample.sce@assays@data$counts
    
    # Subset 
    sample_raw_count <- sample_raw_count[,sample_metadata$spot]
    
    # Rename the count matrix column to the index
    colnames(sample_raw_count) <- sample_metadata$sampleid_spot
    
    sample_count_all <- cbind(sample_count_all, sample_raw_count)
  }
  
  sample_count_all <- sample_count_all[,-1]
  nUMI <- colSums(sample_count_all) # In this case, total counts per pixel is nUMI
  
  # Create SpatialRNA object
  puck <- SpatialRNA(counts=sample_count_all, nUMI=nUMI, use_fake_coords = T)
  message("Original spatial data: ", puck@counts@Dim[1], "x", puck@counts@Dim[2]) 
  
  # QC puck file
  qc_barcodes <- fread(file.path(opt$qc_spots))$sampleid_spot
  qc_genes <- fread(file.path(opt$qc_genes))$gene_name
  puck <- puck %>%
    spacexr::restrict_puck(barcodes = qc_barcodes) 
  
  qc_gene_present <- qc_genes[qc_genes %in% puck@counts@Dimnames[[1]]]
  puck <- spacexr::restrict_counts(puck,
                                   gene_list = qc_gene_present, UMI_thresh = 0, UMI_max = Inf, counts_thresh = 0)
 

  message("QCed spatial data: ", puck@counts@Dim[1], "x", puck@counts@Dim[2])
  
  # Save puck
  saveRDS(puck, file.path(out_dir,paste0("pooled_", 
                                         paste0(region, "_", paste0(unlist(grps), collapse = "_")), 
                                         "_puck.rds")))
  
  
  # Check the distribution of normalized mean gene expression
  count_norm <- sweep(puck@counts, 2, puck@nUMI, "/")
  gene_mean <- Matrix::rowMeans(count_norm)
  
  ggplot(data.frame(gene_mean=gene_mean) %>% rownames_to_column("gene"),
         aes(x=gene_mean))+
    geom_histogram(bins = 1000)+
    scale_x_continuous(trans = "log10")
  scale_x_continuous(limits = c(0, 50))
  
  ggsave(file.path(out_dir,paste0("pooled_", 
                                  paste0(region, "_", paste0(unlist(grps), collapse = "_")), 
                                  "_gene_mean_distribution.png")))  
  
  
  # Create RCTD object
  myRCTD <- create.RCTD(spatialRNA = puck,
                        reference = reference,
                        max_cores = max_core,
                        CELL_MIN_INSTANCE = 0)
  
  
  saveRDS(myRCTD, file.path(out_dir, paste0("rctd_create_pool_", 
                                            paste0(region, "_", paste0(unlist(grps), collapse = "_")), 
                                            "_full_no_filter_c2l.rds")))
  

  rm(reference)
  gc()
  
  # Load customer weight matrix
  c2l_wts <- fread(file.path(opt$ct_wt)) %>%
    column_to_rownames("V1") 

  # # replace NAs in weights with 0
  c2l_wts[is.na(c2l_wts)] <- 0
  
  # Scale the weight matrix to have row sum of 1
  c2l_wts_scaled_mtx <- sweep(c2l_wts, 1, rowSums(c2l_wts), '/')
  
  # Import weight 
  barcodes <- myRCTD@spatialRNA@counts@Dimnames[[2]] 
  c2l_wt_rep <- c2l_wts_scaled_mtx %>% 
    dplyr::filter(rownames(.) %in% barcodes) 
  
  myRCTD <- import_weights(myRCTD, c2l_wt_rep)
  
  # need to set the doublet_mode as "full" so that the code will use this weights in slot: myRCTD@results$weights
  myRCTD@config$doublet_mode <- "full"
  myRCTD@config$RCTDmode <- "full"
  
  # create X with covariate: animal_id
  animal_ID <- st_meta_harmonized %>% 
    dplyr::filter(sampleid_spot %in% barcodes) %>%  
    mutate(sampleid=str_extract(sampleid_spot, ".*?(?=_)")) %>% 
    dplyr::select(type, subject_id, sampleid) %>% 
    mutate(subject_id_n=str_extract(subject_id, ".*?(?=_)")) %>%
    mutate(type=factor(type, levels=c(grps[2],grps[1])))
  
  X <- model.matrix(~ type + type:subject_id_n, data = animal_ID)
  rownames(X) <- barcodes 
  
  
  myRCTD <- run.CSIDE(myRCTD,
                      X,
                      barcodes,
                      cell_types =  cell_types,
                      weight_threshold = weight_thresh,
                      log_fc_thresh = log_fc_thresh,
                      doublet_mode = FALSE, #Ran with full weights, not doublet
                      gene_threshold = gene_threshold,
                      cell_type_threshold = cell_type_count_thresh,
                      # replicate_index = index_tmp, #Set this for CPZ vs CTL replicate groups
                      fdr = 0.05)
  
  
  saveRDS(myRCTD, file.path(out_dir, paste0("cside_pool_", 
                                            paste0(region, "_", paste0(unlist(grps), collapse = "_")), 
                                            "_full_no_filter_c2l.rds")))
  

  res <- myRCTD@de_results$all_gene_list %>% 
    map(function(x) x %>% rownames_to_column("gene")) %>% 
    rbindlist(idcol = "cell") 
  
  write.csv(res, file.path(out_dir, paste0("DEG_pool_", 
                                           paste0(region, "_", paste0(unlist(grps), collapse = "_")), 
                                           "_full_no_filter_c2l.csv")))
}

