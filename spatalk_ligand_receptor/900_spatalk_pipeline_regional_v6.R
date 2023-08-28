#Version history
#v1 basic pipeline with regional selection for section
#v2 added res_out, added manual LR option
#v3 added option for multiple celltype pairs
#v4 using downsample keep count matrix row order and weight matrix column order
#v5 intermediate attempt to put seed before each spatalk call
#v6 modify to include seed in parallel dec_celltype, dec_cci_all, dec_cci

#Load libraries
library(tidyverse)
library(data.table)
library(glue)
library(SpaTalk) #local
library(Seurat)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(

  # ST data input options
  make_option("--st_meta_data", action = "store", default = NA, type = "character",
              help = "the ST meta data [required]"),
  make_option("--qc_spots", action = "store", default = NA, type = "character",
              help = "the QCed spots [required]"),
  make_option("--qc_genes", action = "store", default = NA, type = "character",
              help = "the QCed genes [required]"),
  make_option("--section", action = "store", default = NA, type = "character",
              help = "the section ID to analyze [required]"),
  # section regional subset
  make_option("--region_column", action="store", default=NA, type="character",
              help="spot-level metadata column that describes the grouping"),
  make_option("--region_subset", action="store", default=NA, type="character",
              help="values in region column to subset spots to, separated by space or comma or semicolon."),
  
  # deconvolution weight matrix
  make_option("--ct_wt", action = "store", default = NA, type = "character",
              help = "full path to precomputed ST deconvoluted cell type weight matrix [required]"),
  
  # single cell reference input options
  make_option("--sc_ref_10x_dir", action = "store", default = NA, type = "character",
              help = "the directory of single cell reference file of 10x [required]"),
  make_option("--sc_celltype_meta_data", action = "store", default = NA, type = "character",
              help = "Full path to single cell metadata of cells and cell types given in sc reference --sc_ref_10x_dir [required]"),  
  
  # single cell reference filtering options - downsampling
  make_option("--sc_downsample", action = "store_true", default = FALSE, type = "logical",
              help = "include flag to load down sampling of sc reference (requires spacexr)"),
  make_option("--sc_downsample_path", action = "store", default = NA, type = "character",
              help = "Full path to downsample reference rds file (by default require a spacexr object)."),
  
  # single cell reference filtering options - flat filter
  make_option("--sc_flat_filter", action = "store_true", default = FALSE, type = "logical",
              help = "include flag to run flat filtering of sc reference file"),
  make_option("--gene_min_count", action="store", default=NA, type="numeric",
              help = "minimum count for a gene, used in conjunction with --gene_min_cell"),
  make_option("--gene_min_cell", action="store", default=NA, type="numeric",
              help = "keep gene with minimum number of cells with at least --gene_min_count"),
  make_option("--cell_min_count", action="store", default=NA, type="numeric",
              help = "minimum count per gene for a cell, used in conjunction with --cell_min_gene"),
  make_option("--cell_min_gene", action="store", default=NA, type="numeric",
              help = "keep cell with minimum number of genes with at least --cell_min_count"),
  
  #Step 1 Decon options
  make_option("--ct_decon", action = "store_true", default = FALSE, type = "logical",
              help = "run the deconvolution function [required]"), 
  make_option("--n_perm", action = "store", default = 1000, type = "integer",
              help = "number of permutations [required]"), 
  make_option("--min_percent", action="store", default = 0.5, type="numeric",
              help="minimum proportion of cell type in a spot to retain.  Lower retains less common cell types."),
  
  #Step 2 ligand-receptor analysis options
  make_option("--cci", action = "store_true", default = FALSE, type = "logical",
              help = "run the cell-cell interaction function [required]"), 
  make_option("--if_cci_all", action = "store_true", default = FALSE, type = "logical",
              help = "run the cell-cell interaction for all pairs [required]"),
  make_option("--man_sender_receiver", default = NA, type = "character",
              help = "csv file to run batch of sender and receiver pairs.  Columns are sender, receiver."),
  make_option("--sender", action = "store", default = NA, type = "character",
              help = "sender to specify when running individual sender-receiver pair; required when --cci=TRUE but not --if_cci_all=TRUE or when using --man_sender_receiver [required]"),
  make_option("--receiver", action = "store", default = NA, type = "character",
              help = "sender to specify when running individual sender-receiver pair; required when --cci=TRUE but not --if_cci_all=TRUE or when using --man_sender_receiver [required]"),
  make_option("--man_lr_path", action = "store", default = NA, type = "character",
              help = "include path to csv file to use manual ligand-receptor pairs. Columns are ligand,receptor,species."),
  
  # Step 3 
  make_option("--res_out", action = "store_true", default = FALSE, type = "logical",
              help = "include this step to run the result extraction"),
  make_option("--res_path", action = "store", default = NA, type = "character",
              help = "full path to rds file from step 2 to summarize (could be a cell type pair or all)"),
  
  # output options
  make_option("--out_dir", action = "store", default = NA, type = "character",
              help = "the output directory [required]"),
  
  # computation options
  make_option("--max_core", action = "store", default = 4, type = "integer",
              help = "number of cores in parallel computing [required]. Used in steps 1 and 3."),
  make_option("--seed", action="store", default = 2023, type = "integer",
              help = "seed for analysis used in steps 1 and 2.")
  
)))

# debug
if(FALSE){
  #Directories
  prog_dir="/edgehpc/dept/compbio/projects/SingleCell/programs/dev/spatial_tx/spatalk"
  #ST data directories
  st_dir="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS"
  st_qc_dir="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/443_ST_QC_summary"
  #ST deconvolution directory
  c2l_dir="/edgehpc/dept/compbio/projects/TST11523/2023/cell2loc/result"
  #Single cell reference data directories
  sc_dir="/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed"
  #Manual LR directory
  man_lr_dir="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/spatalk"
  #Output directory
  out_dir="/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/spatalk/901_spatalk_pipeline_regional_v2_cortex"
  
  opt=list(
    st_meta_data = file.path(st_dir,"metadata.tsv"),
    qc_spots = file.path(st_qc_dir,"443_ST_QC_scenario5_spots.csv"),
    qc_genes = file.path(st_qc_dir,"443_ST_QC_scenario5_genes.csv"),
    section = "072A",
    region_column = "annotation",
    region_subset = "cortex",
    ct_wt = file.path(c2l_dir,"cell_type_weight_mtx_all.csv"),
    sc_ref_10x_dir = file.path(sc_dir),
    sc_celltype_meta_data = file.path(sc_dir,"meta_format_clean.csv"),
    sc_downsample = FALSE,
    sc_flat_filter = TRUE,
    gene_min_count = 3,
    gene_min_cell = 10,
    cell_min_count = 5,
    cell_min_gene = 1000,
    ct_decon = TRUE,
    n_perm = 1000,
    min_percent = 0.2,
    cci = TRUE,
    sender = "microglia",
    receiver = "astrocytes",
    man_lr_path = file.path(man_lr_dir,"805_manual_lrs.csv"),
    res_out = TRUE,
    res_path = file.path(out_dir, "072A_spatalk_dec_celltype_cci_microglia_astrocytes.rds"),
    out_dir = file.path(out_dir),
    max_core = 12,
    seed = 2023
           )
}

#Initial setup
message("Setting up section options at:", Sys.time())
out_dir=file.path(opt$out_dir)
section_id <- opt$section

#Set analysis seed
set.seed(opt$seed)
message("Set a seed: ", opt$seed)

#Prevent Rplots.pdf from being generated when running script
#Let window open when interactive
if(!interactive()) {
  pdf(NULL)
}

#Source modified spatalk source code with seeding
source("/edgehpc/dept/compbio/projects/SingleCell/programs/dev/spatial_tx/spatalk/000_spatalk_aux.R")

#Set up an analysis directory with global permissions
if(! dir.exists(out_dir)){
  dir.create(out_dir, recursive=T)
  system(glue("chmod 775 {out_dir}"))
}

#Step 1: Prepare data
if(opt$ct_decon) {
  
  ################# Prepare ST data ###################
  
  #Load the st metadata and join to c2l estimated cell count matrix
  st_meta <- fread(file.path(opt$st_meta_data)) %>%
    mutate(sampleid_spot = paste0(sample,"_",gsub("-.*$","",V1),"-1"),
           bio_rep = paste0(type,"_",BiologicalSampleID),
           cell=sampleid_spot) 
  
  #Load st section names and select section
  sections <- list.files(path="/edgehpc/dept/compbio/users/ezhao/data-raw",full.names=T)[1:28]
  section <- sections[str_detect(sections, section_id)]
  message("Loading section data")
  
  #Load single section
  srt <- Seurat::Load10X_Spatial(section,slice=glue("slice_{section_id}"), use.names=F ) #Do not use the gene names, there are unexplained duplications
  gene_info_raw <- fread(glue("/edgehpc/dept/compbio/users/ezhao/data-raw/{section_id}/filtered_feature_bc_matrix/features.tsv.gz"), header=F)
  fix_names <- scater::uniquifyFeatureNames(ID = rownames(srt), names = gene_info_raw$V2) #Use the same approach as Edward Zhao in BayesSpace
  rownames(srt@assays$Spatial@counts) <- fix_names
  rownames(srt@assays$Spatial@data) <- fix_names
  
  #Load st spot-level QC files
  qc_spots <- fread(file.path(opt$qc_spots))
  qc_genes <- fread(file.path(opt$qc_genes))
  
  # #Load SpaGCN metadata (not used for now)
  # spagcn_map_df <- fread(file.path(qc_dir, "443a_spagcn_17cluster_map.csv"))
  # #Join to existing
  # st_meta <- meta_data # %>% left_join(spagcn_map_df)
  
  #QC ST count matrix
  st_count <- srt@assays$Spatial@counts
  colnames(st_count) <- paste0(section_id,"_",colnames(st_count)) #Need sampleid_spot
  st_count_qc <- st_count[rownames(st_count) %in% qc_genes$gene_name, colnames(st_count) %in% qc_spots$sampleid_spot]
  
  #format region vector
  st_regions <- unlist(strsplit(opt$region_subset, split = "[ ,;]+"))
  
  #subset metadata to region
  st_meta_region <- st_meta %>% filter(!!sym(opt$region_column) %in% st_regions)
  
  #subset spots to region
  st_count_qc <- st_count_qc[,colnames(st_count_qc) %in% st_meta_region$sampleid_spot, drop=FALSE]
  
  #QC and format metadata
  st_meta_qc <- st_meta_region %>% filter(sampleid_spot %in% colnames(st_count_qc)) %>%
    dplyr::select(spot = sampleid_spot,
                  x = array_col,
                  y = array_row)
  
  #Revise the count matrix gene symbols according to the SpaTalk recommendation
  st_count_qc_rev <- rev_gene(data = st_count_qc,
                              data_type = "count",
                              species = "Mouse",
                              geneinfo = geneinfo)
  
  message("Spatial region data after QC: ", st_count_qc_rev@Dim[1], " genes x ", st_count_qc_rev@Dim[2], " spots.")
  
  # create ST SpaTalk data
  obj <- createSpaTalk(st_data = st_count_qc_rev,
                       st_meta = st_meta_qc, #Can only have the three columns
                       species = "Mouse",
                       if_st_is_sc = F, #The Visium data is NOT cell-level
                       spot_max_cell = 30 #Recommended for 10X Visium 55 micron spot
  
  )
  
  ################# Prepare precomputed ST deconvolution result ###################
  
  #Load the c2l estimated cell count matrix assuming V1 contains format V1 column sampleid_UMI
  c2l_df <- fread(file.path(opt$ct_wt)) %>%
    # remove the cell type of "" and not_assigned
    dplyr::select(-c(nan,not_assigned)) %>%
    dplyr::mutate(spot = gsub("\\-","_",V1) ) %>% #Note: SpaTalk does this silently!!!!!
    dplyr::select(-V1) %>%
    dplyr::select(spot, everything() )
  
  #replace NAs in weights with 0
  c2l_df[is.na(c2l_df)] <- 0
  c2l_mat <- as.matrix(c2l_df %>% column_to_rownames("spot"))
  
  #QC the c2l_mat using the ST data
  c2l_mat_qc <- c2l_mat[rownames(c2l_mat) %in% colnames(obj@data$rawdata), ]
  
  ################# Prepare single cell reference ###################
  
  #Load single cell metadata
  sc_meta_data  <-  fread(file.path(opt$sc_celltype_meta_data), data.table=F) %>%
    dplyr::filter(!Celltype_annotation_clean %in% c("x","not_assigned"))
  message("Removed from scRNA-seq reference 7741 cells without celltype and 9380 cells with celltype not assigned from single cell metadata.")
  
  message("Loading full single cell reference count file")
  # load single cell ref
  genes <- fread(file.path(opt$sc_ref_10x_dir, "genes.tsv"), header = F)
  barcodes <- fread(file.path(opt$sc_ref_10x_dir, "barcodes.tsv"), header = F)
  sc_data <- Matrix::readMM(file.path(opt$sc_ref_10x_dir, "matrix.mtx"))
  colnames(sc_data) <- barcodes$V1
  rownames(sc_data) <- genes$V1
  sc_data <- methods::as(sc_data, "dgCMatrix")
  
  #Subset single cell reference matrix to remove unassigned cells
  sc_data <- sc_data[,sc_meta_data$cell]
  
  #Use pre-computed downsampled SC reference
  if(opt$sc_downsample){
    message("using downsampled reference created in spacexr package!  This data will be used instead of full single cell reference.")
    require(spacexr)
    downsample_spacexr <- readRDS(opt$sc_downsample_path)
    sc_data <- downsample_spacexr@counts
    message("Single cell downsampled reference data: ", sc_data@Dim[1], " genes x ", sc_data@Dim[2], " cells.")
  }

  #Perform single-cell flat filtering options
  if(opt$sc_flat_filter){
  message("Performing flat filtering of single cell data")
  message(glue("QC genes: at least {opt$gene_min_count} count in at least {opt$gene_min_cell} cells"))
  # option to filter sc_data reference by genes and cells
  sc_rs <- which(rowSums(sc_data >= opt$gene_min_count) >= opt$gene_min_cell)
  sc_data <- sc_data[sc_rs,]
  message(glue("QC cells: at least {opt$cell_min_count} counts in at least {opt$cell_min_gene} genes"))
  sc_cs <- which(colSums(sc_data >= opt$cell_min_count) >= opt$cell_min_gene)
  sc_data <- sc_data[,sc_cs]
  message("Single cell reference data after QC: ", sc_data@Dim[1], " genes x ", sc_data@Dim[2], " cells.")
  }
  
  #Create vector of cell types
  sc_celltype_df  <-   sc_meta_data %>%
    dplyr::filter(cell %in% sc_data@Dimnames[[2]])
  
  #vector of cell type labels
  sc_celltype <- sc_celltype_df %>%
    pull(Celltype_annotation_clean)

  #Order sc counts same as labels
  sc_data <- sc_data[,sc_celltype_df$cell]
  
  # refresh c2l weight matrix with remaining post-filtered cell types
  # Preserve order of weight matrix columns
  c2l_mat_qc_cols <- colnames(c2l_mat_qc)
  c2l_mat_qc_cols_select <- c2l_mat_qc_cols[which(c2l_mat_qc_cols %in% unique(sc_celltype))]
  
  # refresh c2l weight matrix with remaining post-filtered cell types
  c2l_mat_qc <- c2l_mat_qc[,c2l_mat_qc_cols_select]
  
  ################# Set up SpaTalk object (ST data,  precomputed ST deconvolution result, and single cell reference data) ###################
  
  # Add cell2location 
  # loading the results directly
  # cell type deconvolution using spatalk
  
  # param env When method set to 6, namely use stereoscope python package to deconvolute, please define the python environment of installed stereoscope.
  # Default is the 'base' environment. Anaconda is recommended. When method set to 7, namely use cell2location python package to deconvolute,
  # please install cell2location to "base" environment.
  
  # use_python("/home/spiya/.conda/envs/cell2loc_env/bin/python")
  message("Beginning spatalk ct_decon object step ", Sys.time())
  obj <- SpaTalk::dec_celltype(obj,
                      sc_data = sc_data, #The error handling for this function requires this to be defined
                      sc_celltype = sc_celltype,
                      iter_num = opt$n_perm, # default 1000
                      use_n_cores = opt$max_core,
                      min_percent = opt$min_percent,
                      # method = 7, # 1 for spatalk, 2 for RCTD, 3 for Seurat, 7 for cell2location
                      # env = "/home/spiya/.conda/envs/cell2loc_env",
                      dec_result = c2l_mat_qc, # will not perform deconvolution
                      rand_seed = opt$seed) 
  
  saveRDS(obj, file.path(out_dir, paste0(section_id, "_c2l_wt_mtx", glue("_spatalk_dec_celltype.rds") )))
  
  message("Finished ct_decon step at ", Sys.time())
  
}

#Step 2: SpaTalk analysis 
if(opt$cci) {
  obj <- readRDS(file.path(out_dir, paste0(section_id, "_c2l_wt_mtx","_spatalk_dec_celltype.rds")))
  
  message("Start find_lr_path analysis at: ",Sys.time())
  data("lrpairs")
  data("pathways")
  
  if(!is.na(opt$man_lr_path)){
    lr_add_df <- fread(opt$man_lr_path, data.table=F)
  } else {
    lr_add_df <- NULL
  }
  
  #Bind the paths to the SpaTalk data
  lrpairs <- rbind(lrpairs,lr_add_df) %>% distinct()
  
  obj <- find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
  
  if( opt$if_cci_all ){
    message("Start cci analysis at: ",Sys.time())
    message("Analyzing all sender-receiver pairs...")
    obj <- SpaTalk::dec_cci_all(obj,
                       n_neighbor = 10, # default
                       min_pairs = 5, # default
                       min_pairs_ratio = 0, # default
                       per_num = 1000, # default
                       co_exp_ratio = 0.1, # default
                       use_n_cores = opt$max_core, 
                       rand_seed = opt$seed #Provide seed for parallel
                       )
    
    saveRDS(obj, file.path(out_dir, paste0(section_id, 
                                           "_c2l_wt_mtx",
                                           "_spatalk_dec_celltype_cci_all.rds")))
    
    message("Finished all cci at ", Sys.time()) 
    
  } else {
    
    #Either run a batch or run a single pair
    if(!is.na(opt$man_sender_receiver)){
      man_sender_receiver <- fread(opt$man_sender_receiver, data.table=F)
      single_sender_receiver <- NULL
    } else{
      #Flags for sender and receiver must be defined if not using --if_cci_all or --man_sender_receiver
      stopifnot(exprs = {!is.na(opt$sender)
        !is.na(opt$receiver)} )
      man_sender_receiver <- NULL
      single_sender_receiver <- data.frame(sender=opt$sender, receiver=opt$receiver)
    }
    
    #Multiple pairs are allowed in the data.frame
    sender_receiver_df <- rbind(man_sender_receiver, single_sender_receiver) %>% distinct()
    
    #Run the pairings
    #Error if pairing does not work
    for(i in 1:nrow(sender_receiver_df)){
    message("Start cci analysis at: ",Sys.time())
    message("Analyzing sender-receiver pair: ", sender_receiver_df$sender[[i]], " ", sender_receiver_df$receiver[[i]])
    
    tryCatch({
    obj <- SpaTalk::dec_cci(obj,
                   celltype_sender = sender_receiver_df$sender[[i]],
                   celltype_receiver = sender_receiver_df$receiver[[i]],
                   n_neighbor = 10, # default
                   pvalue = 0.05, # default
                   min_pairs = 5, # default 5
                   min_pairs_ratio = 0, # default
                   per_num = 1000, # default 1000
                   co_exp_ratio = 0.1, # default
                   if_doParallel = TRUE, # default  
                   use_n_cores = opt$max_core,
                   rand_seed = opt$seed)}, 
    error=function(e, obj_err=obj){
      message("Something went wrong for this pairing.")
      print(e)
      return(obj_err)
      }
    
    )
    
    }
    
    #Save this file out as an rds 
    saveRDS(obj, file.path(out_dir, paste0(section_id, 
                                           "_c2l_wt_mtx",
                                           glue("_spatalk_dec_celltype_cci.rds"))))
    
    message("Finished cci at ", Sys.time()) 
    
  }
  
}

#Step 3: Extract analysis results
if(opt$res_out){
  
  message("started --res_out at ", Sys.time() )
  
  #Read in step 2 analysis object
  obj <- readRDS(file.path(opt$res_path))
  
  #Set up output name
  res_out_name <- gsub("\\.rds$","",basename(opt$res_path))
  
  #Extract LR interaction result table
  lri <- obj@lrpair
  
  #Write out step 2 LR permutation test results
  data.table::fwrite(lri, file.path(out_dir, glue("{res_out_name}_lri.csv") ) )
  
  #Prepare LR pathway result
  #Throw away distance matrix
  obj2 <- obj
  obj2@dist <- matrix()
  #Run in parallel the rows of the LR permutation test results
  message("Performing pathway analysis (warning: can be quite slow for lots of cell types)")
  lr_path_list <- mclapply(1:nrow(lri), function(x, obj_input = obj2){
    
    #Isolate significant LRs
    lri_input <- obj_input@lrpair
    
    #Perform pathway analysis
    res <- get_lr_path(object = obj_input,
                celltype_sender = lri_input$celltype_sender[[x]],
                celltype_receiver = lri_input$celltype_receiver[[x]],
                ligand = lri_input$ligand[[x]],
                receptor = lri_input$receptor[[x]],
                min_gene_num = 3 
                )
    
    return(res)
    
  },mc.cores = opt$max_core)
  
  #Bind the tf_path results
  tf_df <- rbindlist( lapply(lr_path_list, "[[", 1) ) %>% distinct()
  
  #Bind the path_value results
  lr_path_df <- rbindlist( lapply(lr_path_list, "[[", 2) )
  
  #Write out tf results 
  write.csv(tf_df, file.path(out_dir, glue("{res_out_name}_tf_distinct.csv") ) )
  
  #Write out pathway results
  write.csv(lr_path_df, file.path(out_dir, glue("{res_out_name}_path.csv") ) )
  
  message("Finished basic outputs at ", Sys.time() )
  
  #Graphical outputs depend on the sender and receiver pairs measured in the outputs above
  
  #Create dataframe of sender receiver cell types that were measured in results
  lri_celltype_df <- lri %>%
    dplyr::select(sender=celltype_sender,receiver=celltype_receiver) %>% distinct()
    
  #Output graphical results
  mclapply(1:nrow(lri_celltype_df), function(x, obj_input=obj, lri_celltype_df_input=lri_celltype_df){
  
  tryCatch({
  #Plot the made-up spatial relationships that spatalk created
  ccdist_plot <- plot_ccdist(object = obj_input, 
              celltype_sender = lri_celltype_df_input$sender[[x]], 
              celltype_receiver = lri_celltype_df_input$receiver[[x]],            
              size = 0.2,
              arrow_length = 0.02)
  
  ggsave(plot = ccdist_plot, 
         filename = file.path(out_dir, glue("{section_id}_{lri_celltype_df_input$sender[[x]]}_{lri_celltype_df_input$receiver[[x]]}_ccdist.png") ), 
         width=8,height=6,units="in",dpi=300)},
  error=function(e){
    message("Could not write out distance plot")
    print(e)
  })
  
  #If a top 20 plot exists, print it
  tryCatch({
  
  top20_lri <- plot_cci_lrpairs(object = obj_input,
                                celltype_sender = lri_celltype_df_input$sender[[x]],
                                celltype_receiver = lri_celltype_df_input$receiver[[x]],
                                type = "sig")
  
  ggsave(plot = top20_lri, 
         filename = file.path(out_dir, glue("{section_id}_{lri_celltype_df_input$sender[[x]]}_{lri_celltype_df_input$receiver[[x]]}_top20_cci_lrpairs.png") ), 
         width=5,height=5,units="in",dpi=300)
  
  }, error=function(e){
    message("Could not write out top 20 plot.")
    print(e)
  })
  
  }, mc.cores = opt$max_core)
  
  message("Finished outputting graphical results at ", Sys.time() )
  
}
