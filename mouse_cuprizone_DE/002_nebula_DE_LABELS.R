library(tidyverse)
library(data.table)
library(glue)
library(Matrix)
library(optparse)
library(ggrepel)
library(clusterProfiler)

# User must set up a program directory
prg.dir <- file.path(getwd()) #"/mnt/depts/dept04/compbio/projects/SingleCell/programs/dev/mouse_cuprizone_DE_interactivereport"

# Make options
option_list = list(
  make_option("--in_dir", action = "store", default = NA, type = "character",
              help = "Path to input files [required]"),
  make_option("--gene_info_file", action = "store", default = NA, type = "character",
              help = "File name of the gene info file [required]"),
  make_option("--meta_file", action = "store", default = NA, type = "character",
              help = "File name of the meta file [required]"),
  make_option("--count_file", action = "store", default = NA, type = "character",
              help = "File name of the count data [required]"),
  make_option("--celltype_col",action="store",default=NA,type="character",
              help="celltype column in metadata"),
  make_option("--sample_col",action="store",default=NA,type="character",
              help="sample column in metadata"),
  make_option("--treatment_col",action="store",default=NA,type="character",
              help="treatment column in metadata"),
  make_option("--cluster", action = "store", default = "cell_type", type = "character",
              help = "Which cluster to use? [default: %default]"),
  make_option("--covars", action = "store", default = NA, type = "character",
              help = "The names of covars to be included for the analysis"),
  make_option("--ref_group",action="store",default=NA,type="character",
              help="treatment contrast reference level."),
  make_option("--alt_group",action="store",default=NA,type="character",
              help="treatment contrast alternative level."),
  make_option("--lib_size_high", action="store", default=20E6, type="numeric",
              help = "maximum library size per cell for 1st round filtering in R6 pipeline."),
  make_option("--minimum.cells.per.gene.type", action = "store", default = "and", type = "character",
              help = "and: both groups must have at least 50 cells, or: either group (or both) must have at least 50 cells [default: %default]"),
  make_option("--minimum.cells.per.subject", action = "store", default = 5, type = "integer",
              help = "Minimum number of cells per subject [default: %default]"),
  make_option("--out_dir", action = "store", default = NA, type = "character",
              help = "Path for output file [required]"),
  make_option("--out_prefix", action = "store", default = NA, type = "character",
              help = "Prefix of the output files [required]"),
  make_option("--volcano", action="store_true", default=FALSE, type="logical",
              help = "Set TRUE to output volcano plot of DE result.")
)

# Configure optparse 
opt = parse_args(OptionParser(option_list=option_list))
print(opt)

if (is.na(opt$covars)) {
  covars = NULL
} else {
  covars = unlist(strsplit(opt$covars, split = "[ ,;]+"))
}

if(file.exists( file.path(opt$out_dir, glue("{opt$out_prefix}_{opt$cluster}.csv") ) ) ){
  message("DE analysis already performed, reading local file.")
  Nebula.HL_results <- read.csv( file.path(opt$out_dir, glue("{opt$out_prefix}_{opt$cluster}.csv") ) )
} else {
# Read in data
UMI_data<-readMM(file.path(opt$in_dir, opt$count_file)) # assumes a dgTMatrix in the correct shape
meta_info=read.csv(file.path(opt$meta_file), stringsAsFactors=FALSE)
gene_info=fread(file.path(opt$in_dir, opt$gene_info_file), stringsAsFactors=FALSE, header=F) %>% #saved out as a headerless tsv w/2 columns
  data.frame()
  
colnames(UMI_data) <- meta_info$cell
rownames(UMI_data) <- gene_info[,1] 

# Load the pipeline code
source(file.path(prg.dir, "pipeline_class.R"))

#Load file
sce <- BiostatsSingleCell$new(count_data = UMI_data,
                              meta_data = meta_info,
                              sampleId_col = opt$sample_col,
                              cluster_col = opt$celltype_col,
                              treatment_col = opt$treatment_col)

# Filtering round 1
sce$apply_filter(min.perc.cells.per.gene = 0.00, lib_size_high = opt$lib_size_high)

# Apply group mode
sce$set_group_mode(cluster_of_interest = opt$cluster, ref_group = opt$ref_group, alt_group = opt$alt_group)

# Filtering round 2
sce_qc <- sce$apply_filter_contrasts_R6(min.cells.per.gene.type=opt$minimum.cells.per.gene.type, 
                                        min.perc.cells.per.gene=0.1, 
                                        min.cells.per.subj=opt$minimum.cells.per.subject)

#Run NEBULA HL
Nebula.HL_results = sce_qc$nebula_pipeline(covs = covars,method="HL")$res.tab

#Write out results
write.csv(Nebula.HL_results, file.path(opt$out_dir, glue("{opt$out_prefix}_{opt$cluster}.csv")), row.names = FALSE)
}

if(opt$volcano){
  
  source(file.path(prg.dir,"000_nebula_DE_func","volcano_func.R"))
  
  volc_plot <- volcano_func(result_df = Nebula.HL_results,FDR_threshold = 0.05,FC_threshold = 1.5)
  
  ggsave(plot = volc_plot, 
         filename = file.path(opt$out_dir, glue("{opt$out_prefix}_{opt$cluster}_volcano.png")), 
         width=7.5, height=5, units="in", dpi=300)
  
}
