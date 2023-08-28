require(tidyverse)
require(magrittr)
library(data.table)
library(Matrix)
library(janitor)

#Set the directory
data_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_ambient_doublet_removed")
out_dir <- data_dir

#Load
meta_data_format <- fread(file.path(out_dir, "meta_format.tsv"))

#Clean cell type names as in the NEBULA DE 
#Do some basic metadata cleaning
#Fix the cell type labels
clean_annotations <- data.frame(Celltype_annotation = unique(meta_data_format$Celltype_annotation),
                                Celltype_annotation_clean = janitor::make_clean_names(unique(meta_data_format$Celltype_annotation)), check.names = F 
)
#Add clean cell type labels to metadata
meta_data_format_clean <- meta_data_format %>%
  left_join(clean_annotations) %>% data.frame

fwrite(meta_data_format_clean, file.path(out_dir, "meta_format_clean.csv"))
