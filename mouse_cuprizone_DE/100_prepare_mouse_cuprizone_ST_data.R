library(tidyverse)
library(data.table)
library(glue)
library(Matrix)
library(optparse)
library(ggrepel)
library(clusterProfiler)
library(Seurat)
library(SingleCellExperiment)

#Set up directories
data_dir <- "/edgehpc/dept/compbio/users/ezhao/out/harmony"
out_dir <- "/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE/CRZ_vs_CTRL_ST/"

#Read in QC'ed rds data
#The thresholds used here were spot-level per sample, removed if < 100 reads or > 45% mitochondrial reads
#This data is read in as a SingleCellExperiment file
st_data <- readRDS(file = file.path(data_dir, "sce.preprocessed.RDS"))

#Extract metadata and format the index
metadata_qc <- data.frame(colData(st_data),stringsAsFactors = FALSE) %>%
  mutate(sampleid_spot = paste0(SampleID,"_",spot)) #This value must be unique

library(BayesSpace, lib.loc="/home/ezhao1/R/x86_64-pc-linux-gnu-library/4.1")
#Load data
#Reads data from Space Ranger output into R as a SingleCellExperiment
set.seed(100)
ID.map = read.csv("/edgehpc/dept/compbio/users/ezhao/data-raw/ID.map.csv")
ID.map$slide = substr(ID.map$SampleID, 1,3) #get slide number
ID.map$slide.letter = substr(ID.map$SampleID, 4,4)
method = "harmony"

sample.list = list() #read in all samples
for (i in 1:nrow(ID.map)){
  sample.id = ID.map$SampleID[i]
  sample.list[i] = readVisium(paste0("/edgehpc/dept/compbio/users/ezhao/data-raw/", sample.id))
  sample.list[[i]]$SampleID = sample.id
}
names(sample.list) = ID.map$SampleID

#Combine raw visium data
sce.combined = do.call(cbind, sample.list)

#Extract raw counts and metadata, format the index
st_raw_count <- sce.combined@assays@data$counts
st_raw_metadata <- data.frame(colData(sce.combined),stringsAsFactors = FALSE) %>%
  mutate( sampleid_spot = paste0(SampleID,"_",spot) )

#Ensure correct order of the columns in the raw counts
all.equal(colnames(st_raw_count), st_raw_metadata$spot)
#[1] TRUE

#Rename the count matrix column to the index
colnames(st_raw_count) <- st_raw_metadata$sampleid_spot

#Subset raw metadata to qc sampleid_spots and add in the qc columns
st_qc_metadata <- st_raw_metadata %>% 
  dplyr::filter( sampleid_spot %in% metadata_qc$sampleid_spot) %>%
  left_join(metadata_qc) %>%
  mutate(annotation_index = paste0(spot,"-",SampleID)) #Format difference

#Subset the count matrix to qc sampleid_spots 
st_qc_count <- st_raw_count[, st_qc_metadata$sampleid_spot]

#Ensure correct order of the columns in the raw counts
all.equal(colnames(st_qc_count), st_qc_metadata$sampleid_spot)
#[1] TRUE

#Clean up
rm(sample.list)
rm(st_data)
rm(sce.combined)
rm(st_raw_count)
rm(st_raw_metadata)
rm(metadata_qc)

#Add the metadata for region
st_annotation_meta <- fread(file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/metadata_ST_region.tsv")) %>%
  dplyr::select(annotation_index=V1, annotation)

#Metadata format for seurat
st_qc_metadata_format <- st_qc_metadata %>% 
  left_join(st_annotation_meta) %>%
  column_to_rownames("sampleid_spot") %>%
  mutate(dummy_clust = "dummy")

#Hold the pseudo-single-cell count and metadata as a Seurat object
st_srt_dummy <- Seurat::CreateSeuratObject(counts = st_qc_count, project = "mouse_cuprizone_ST", meta.data = st_qc_metadata_format)

#Create regional srt objects
st_srt_dummy_CC <- subset(st_srt_dummy, subset = annotation == "CC")
st_srt_dummy_cortex <- subset(st_srt_dummy, subset = annotation == "cortex")
st_srt_dummy_hippocampus <- subset(st_srt_dummy, subset = annotation == "hippocampus")
st_srt_dummy_hypoTH <- subset(st_srt_dummy, subset = annotation == "hypoTH")
st_srt_dummy_TH <- subset(st_srt_dummy, subset = annotation == "TH")

#Keep the unassigned spots in a separate seurat object as a backup
st_srt_dummy_unassigned <- subset(st_srt_dummy, subset = annotation == "unassigned")

#Write out the files
saveRDS(st_srt_dummy, file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/mouse_cuprizone_combined_qc_seurat_dummy.rds"))
saveRDS(st_srt_dummy_CC, file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/mouse_cuprizone_combined_qc_seurat_dummy_CC.rds"))
saveRDS(st_srt_dummy_cortex, file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/mouse_cuprizone_combined_qc_seurat_dummy_cortex.rds"))
saveRDS(st_srt_dummy_hippocampus, file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/mouse_cuprizone_combined_qc_seurat_dummy_hippocampus.rds"))
saveRDS(st_srt_dummy_hypoTH, file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/mouse_cuprizone_combined_qc_seurat_dummy_hypoTH.rds"))
saveRDS(st_srt_dummy_TH, file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/mouse_cuprizone_combined_qc_seurat_dummy_TH.rds"))
saveRDS(st_srt_dummy_unassigned, file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/mouse_cuprizone_combined_qc_seurat_dummy_unassigned.rds"))
