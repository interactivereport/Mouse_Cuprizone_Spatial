
prog_dir <- getwd() #Set for code repo
source(file.path("lr_plot.R"))

#load the data
# specify ER and LR
# specify cell_type_pair of interest (Astro and neurons)
basedir <- "/mnt/depts/dept04/compbio/projects/SingleCell/results/dev/400_lr_pipeline/mouse_cuprizone_region_pool_LABELS_truncMean_popsize"
contra <- "CuprizonevsControl"
cellchat <- readRDS(file.path(basedir, contra, "cellchat.rds"))
cell_type_pair2 <- cellchat@meta$final_label %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(cell_type = ".") %>% 
  filter(str_detect(cell_type, "^Ext|^Inh")) %>% 
  mutate(cell_type_pair = paste0(cell_type," -> Micro")) %>% 
  pull(cell_type_pair)

rm(basedir, contra, cellchat)

cell_type_pair1 <- c("Astro -> Micro",
                    "Micro -> Micro")
lr_pair1 <- c("Apoe-Trem2", "Tnf-Tnfrsf1a", "Tnf-Tnfrsf1b")
lr_pair2 <- c("Cx3cl1-Cx3cr1")

lr_pair_list <- list(lr_pair1, lr_pair2)
cell_type_pair_list <- list(cell_type_pair1, cell_type_pair2)

#Run the function on the LR pairs from Hui-Hsin
for (i in 1:2) {
  lr_pair=lr_pair_list[[i]]
  cell_type_pair=cell_type_pair_list[[i]]
  
  # loop for different contrasts
  # set new output file names before each run to avoid overwriting previous results
  wrap_fun <- function(contrast){
    contra <- contrast
    basedir <- "/mnt/depts/dept04/compbio/projects/SingleCell/results/dev/400_lr_pipeline/mouse_cuprizone_region_pool_LABELS_truncMean_popsize"
    out_dir <- file.path(basedir, contra)
    cellchat <- readRDS(file.path(basedir, contra, "cellchat.rds"))
    cell_interest="Micro"
    
    res <- plot_lr(cellchat = cellchat,                    # specify the cellchat object
                   contra = contra,          # specify the contrast
                   cell_interest = cell_interest,          # specify a cell type 
                   source_target = "target",               # determine the cell type as source (emitter) or target (receiver)               
                   LR_pair = lr_pair,
                   cell_type_pair = cell_type_pair,        # emitter or receiver should be the cell_interest parameter
                   pval_cut = 3,                           # plot interactions with pval < 0.01 only
                   regulated = F,                          # need to set the direct_size if TRUE
                   direct_size = 2)                        # size of the direction label
    
    rank_plot_data <- res[["rank_data"]]
    fwrite(rank_plot_data, file.path(out_dir, paste0(cell_interest, "_run_", i, "_emitter_compare_interaction_prob_rank.csv")))
    common_interaction_plot_data <- res[["common_data"]]
    fwrite(common_interaction_plot_data, file.path(out_dir, paste0(cell_interest, "_run_", i, "_LR_compare_common_interaction_receiver.csv")))
    
    p1 <- res[["rank_plot"]]+
      theme(plot.margin=unit(c(1,1,1,1),"cm"),
            plot.background = element_rect(fill="white", linetype = NULL, color = "white")
            )
    ggsave(filename=file.path(out_dir, paste0(cell_interest, "_run_", i, "_emitter_compare_interaction_prob_rank.png")),p1,
           device="png", width=12, height=6, dpi = 300)
    
    p2 <- res[["common_plot"]]+
      theme(plot.margin=unit(c(1,1,1,1),"cm"),
            plot.background = element_rect(fill="white", linetype = NULL, color = "white")
            )
    ggsave(filename=file.path(out_dir, paste0(cell_interest, "_run_", i, "_LR_compare_common_interaction_receiver.png")),p2,
           device="png", width=12, height=10, dpi=300)
    
  }

  map(list("CuprizonevsControl"), wrap_fun)
  
}

