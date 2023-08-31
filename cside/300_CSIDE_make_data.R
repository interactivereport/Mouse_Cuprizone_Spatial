library(tidyverse)
library(data.table)
library(spacexr, lib.loc="/home/jluo/R/x86_64-pc-linux-gnu-library/4.2")
library(doParallel)
library(ggplot2)
library(Seurat)
library(glue)
library(ggpubr)
library(ggrepel)
library(BayesSpace, lib.loc="/home/ezhao1/R/x86_64-pc-linux-gnu-library/4.1")

#Prepare CSIDE formatted data
set.seed(100)
ID.map = read.csv("/edgehpc/dept/compbio/users/ezhao/data-raw/ID.map.csv")
ID.map$slide = substr(ID.map$SampleID, 1,3) #get slide number
ID.map$slide.letter = substr(ID.map$SampleID, 4,4)

#Read harmonized spot data
st_meta_harmonized <- read.csv("/mnt/depts/dept04/compbio/projects/SingleCell/results/dev/spatial_tx/cside/st_meta_harmonized.csv")

#Directory to store puck files for CSIDE annotation contrast analysis
puck_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS/spacexr_format")
system(glue("mkdir -p {puck_dir}"))

#Single region vector
single_annot <- unique(st_meta_harmonized$annotation)
#Remove unassigned
single_annot <- single_annot[which(!grepl("unassigned",single_annot))]

#Pairwise regional vector
pairwise_annot <- apply(combn(unique(st_meta_harmonized$annotation),2),2,paste,collapse='_')
#Remove pairs with unassigned
pairwise_annot <- pairwise_annot[which(!grepl("unassigned",pairwise_annot))]

#Function to make puck files
make_annot_puck <- function(metadata_df, sample_id, annot, out_dir){
  
  #Load using basic input
  sample.sce = readVisium(paste0("/edgehpc/dept/compbio/users/ezhao/data-raw/", sample_id))
  
  #Get the spots from the data to use for each puck
  spot_use <- st_meta_harmonized %>%
    dplyr::filter(annotation %in% annot & grepl(glue("^{sample_id}\\_"), sampleid_spot)) %>%
    pull(sampleid_spot)
  
  #Extract raw counts and metadata, format the index
  if( length(spot_use) > 0 ){
  sample_metadata <- data.frame(colData(sample.sce),stringsAsFactors = FALSE) %>%
    mutate( sampleid_spot = paste0(sample_id,"_",spot) ) %>%
    dplyr::filter(sampleid_spot %in% spot_use)
  
  #Get raw spot count matrix
  sample_raw_count <- sample.sce@assays@data$counts
  
  #subset spots to used spots
  sample_raw_count <- sample_raw_count[,sample_metadata$spot]
  
  #Rename the count matrix column to the index
  colnames(sample_raw_count) <- sample_metadata$sampleid_spot
  
  #Isolate coords
  coords <- sample_metadata %>% dplyr::select(barcodes=sampleid_spot, xcoord=imagecol, ycoord=imagerow) %>%
    rownames_to_column("spot_tmp") %>%
    column_to_rownames("barcodes") %>%
    dplyr::select(-spot_tmp)
  
  #Isolate UMIs
  nUMI <- colSums(sample_raw_count) # In this case, total counts per pixel is nUMI
  
  #Drop spot with 0 nUMI
  nUMI_zero <- nUMI[which(nUMI == 0)]
  message(glue("dropped {length(nUMI_zero)} spots in {sample_id} {annot} with 0 library size"))
  
  nUMI_nz <- nUMI[which(nUMI != 0)]
  coords_nz <- coords %>% 
    rownames_to_column("sampleid_spot") %>% 
    filter(sampleid_spot %in% names(nUMI_nz)) %>% 
    column_to_rownames("sampleid_spot")
  sample_raw_count_nz <- sample_raw_count[,colnames(sample_raw_count) %in% names(nUMI_nz)]
  
  #Create SpatialRNA object
  puck <- SpatialRNA(coords_nz, sample_raw_count_nz, nUMI_nz)
  
  #Write out SpatialRNA object
  saveRDS(puck, file.path(out_dir,glue("{sample_id}_{paste(annot,collapse='_')}_puck.rds") ))
  
  message(glue("Finished {sample_id}_{annot}"))
  } else{
    message(glue("No spots in {sample_id}_{annot}"))
  }
  
}

#Create the grid
puck_df <- setNames(expand.grid(ID.map$SampleID, c(single_annot,pairwise_annot),stringsAsFactors = F), c("SampleID","annot")) %>%
  mutate(SampleID_annot = paste0(SampleID,"_",annot))

#Debug 
if(FALSE){
  
  x <- puck_df$SampleID_annot[250]
  
  #Get sample id
  sample_id <- str_split_fixed(x, pattern="_", n=2)[,1]
  #Get annot
  annot <- str_split_fixed(x, pattern="_", n=2)[,2]
  #Split annot if underscore, otherwise annot
  if(grepl("_",annot)){
    annot <- as.vector(str_split_fixed(annot, pattern="_", n=2))
  }
  
  test <- make_annot_puck(metadata_df = st_meta_harmonized, 
                          sample_id = sample_id,
                          annot = annot,
                          out_dir = puck_dir
                          )
}

#Run rowwise on grid to make the puck files
mclapply(puck_df$SampleID_annot, function(x){
  tryCatch({
  #Get sample id
  sample_id <- str_split_fixed(x, pattern="_", n=2)[,1]
  #Get annot
  annot <- str_split_fixed(x, pattern="_", n=2)[,2]
  #Split annot if underscore, otherwise annot
  if(grepl("_",annot)){
    annot <- as.vector(str_split_fixed(annot, pattern="_", n=2))
  }
  
  make_annot_puck(metadata_df = st_meta_harmonized,
                  sample_id = sample_id,
                  annot = annot,
                  out_dir = puck_dir
                  )
  }, error=function(e) message(glue("Function failed for {sample_id}_{annot}") ) )
}, mc.cores = 4
)

#saveRDS(quiet_list, file.path(puck_dir, "quietly_log.rds"))

#Remove pairwise annot where one of the region is zero spots (there are a few of them)

rm_df <- fread("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS/metadata.tsv") %>% 
  dplyr::count(sample,annotation) %>%
  complete(sample, annotation, fill = list(n = 0)) %>%
  dplyr::filter(n == 0) %>%
  ungroup() %>% group_by(sample) %>%
  summarize(annotation_collapse=paste(annotation, collapse="|")) %>%
  mutate(rm_regex = glue("{sample}.*({annotation_collapse})"))

rm_files <- lapply(rm_df$rm_regex, function(x) list.files(path=puck_dir, pattern=x, full.names = T))

rm_cmd <- lapply(rm_files, function(x) system(glue("rm {paste(x,collapse=' ')}")))
