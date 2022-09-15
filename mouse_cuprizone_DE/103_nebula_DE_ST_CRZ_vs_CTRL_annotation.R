library(tidyverse)
library(data.table)
library(glue)
library(Matrix)
library(optparse)
library(ggrepel)
library(clusterProfiler)
library(Seurat)

# Make options
option_list = list(
  make_option("--annotation", action = "store", default = NA, type = "character",
              help = "annotation region")
)

# Configure optparse 
opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# User must set up paths 
prg.dir <- file.path(getwd()) #"/mnt/depts/dept04/compbio/projects/SingleCell/programs/dev/mouse_cuprizone_DE_interactivereport"
# Load the pipeline code
source(file.path(prg.dir, "pipeline_class.R"))

# Get path
rds_path <- file.path(
  glue("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st/mouse_cuprizone_combined_qc_seurat_dummy_{opt$annotation}.rds")
  )

#Create object from Seurat
sce <- BiostatsSingleCell$new(rds_file = rds_path,
                              sampleId_col = "Number",
                              cluster_col = "dummy_clust", #dummy value
                              treatment_col = "Type" #Uses CPZ/RCV/CTL 
)

#Set up contrast
sce$set_group_mode(cluster_of_interest = "dummy", ref_group = "CTL", alt_group = "CPZ")

#Run NEBULA HL
Nebula.HL_results = sce$nebula_pipeline(covs = NULL,method="HL")$res.tab

#Write out csv
out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE/CRZ_vs_CTRL_ST")
write.csv(Nebula.HL_results, file.path(out_dir, glue("nebula_hl_CRZ_vs_CTRL_ST_{opt$annotation}.csv")), row.names = FALSE)

#Volcano plot
source(file.path(prg.dir,"000_nebula_DE_func","volcano_func.R"))
volc_plot <- volcano_func(result_df = Nebula.HL_results,FDR_threshold = 0.05,FC_threshold = 1.5)
ggsave(plot = volc_plot, 
       filename = file.path(out_dir, glue("nebula_hl_CRZ_vs_CTRL_ST_volcano_{opt$annotation}.png")), 
       width=7.5, height=5, units="in", dpi=300)
