library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratObject)
library(glue)
library(Matrix)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
#Load the metadata with c2l estimated cell count matrix
meta_data <- fread(file.path("/edgehpc/dept/compbio/projects/SingleCell/data/mouse_cuprizone_st_LABELS/metadata.tsv")) %>%
  mutate(sampleid_spot = paste0(sample,"_",gsub("-.*$","",V1),"-1"),
         bio_rep = paste0(type,"_",BiologicalSampleID),
         dummy_clust = "dummy", cell=sampleid_spot)

#Load the raw data 
sections <- list.files(path="/edgehpc/dept/compbio/users/ezhao/data-raw/",full.names=T)[1:28]
section_names <- lapply(sections, basename)
srt_list <- list()

for(i in 1:length(sections)){
  srt_list[[i]] <- Seurat::Load10X_Spatial(sections[i],slice=glue("slice_{unlist(section_names[i])}"), use.names=F ) #Do not use the gene names, there are unexplained duplications
  gene_info_raw <- fread(glue("/edgehpc/dept/compbio/users/ezhao/data-raw/{unlist(section_names[i])}/filtered_feature_bc_matrix/features.tsv.gz"), header=F)
  fix_names <- scater::uniquifyFeatureNames(ID = rownames(srt_list[[i]]), names = gene_info_raw$V2) #Use the same approach as Edward Zhao in BayesSpace 
  rownames(srt_list[[i]]@assays$Spatial@counts) <- fix_names
  rownames(srt_list[[i]]@assays$Spatial@data) <- fix_names
}

rm(gene_info_raw); rm(fix_names)

#Seurat QC steps
qc_func <- function(merged_srt){
  #Calculate features
  merged_srt <- PercentageFeatureSet(merged_srt, "^mt-", col.name="percent_mito")
  merged_srt <- PercentageFeatureSet(merged_srt, "^Hb.*-", col.name="percent_hb")
  merged_srt <- PercentageFeatureSet(merged_srt, "(^Rpl|^Rps)", col.name="percent_ribo")
  #Spot-level filter the same as used for the region-level comparison
  merged_srt_qc = merged_srt[, merged_srt$nCount_Spatial > 100 & merged_srt$percent_mito < 45 ]
  #We can remove mitochondrial genes and Bc1
  merged_srt_qc <- merged_srt_qc[!grepl("(Bc1|^mt-)", rownames(merged_srt_qc)), ]
  #Normalize remaining data
  merged_srt_qc <- SCTransform(merged_srt_qc, assay = "Spatial", verbose = TRUE, method = "poisson")
  return(merged_srt_qc)
}

srt_ctl_cpz_rcv <- merge(x= srt_list[[2]], #072B
                 y= c(
                      #CTL
                      srt_list[[10]],
                      srt_list[[13]],
                      srt_list[[16]],
                      srt_list[[18]],
                      srt_list[[19]],
                      srt_list[[21]],
                      srt_list[[23]],
                      srt_list[[27]],
                      srt_list[[28]],
                      #CPZ
                      srt_list[[1]], #072A
                      srt_list[[4]],
                      srt_list[[7]],
                      srt_list[[8]],
                      srt_list[[12]],
                      srt_list[[15]],
                      srt_list[[17]],
                      srt_list[[22]],
                      srt_list[[25]],
                      #RCV
                      srt_list[[3]], #072C
                      srt_list[[5]],
                      srt_list[[6]],
                      srt_list[[9]],
                      srt_list[[11]],
                      srt_list[[14]],
                      srt_list[[20]],
                      srt_list[[24]],
                      srt_list[[26]]
                      ),
                 #Name the barcodes uniquely
                 add.cell.id =
                   c(
                     #CTL
                     "072B","075B","076A","076D","106B","106C","107A","107C","108C","108D",
                     #CPZ
                     "072A","072D","074C","074D","075D","076C","106A","107B","108A",
                     #RCV
                     "072C","074A","074B","075A","075C","076B","106D","107D","108B"
                                                                   )
                 )

#Clean up the environment
rm(srt_list)
gc()

#Set an out directory
out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/cside/443_ST_QC_summary")

#Load SpaGCN cluster map
spagcn_map <- fread(file.path(out_dir, "443a_spagcn_17cluster_map.csv")) %>% select(cell=sampleid_spot, spagcn_clust) %>%
  mutate(spagcn_clust = as.character(spagcn_clust))

#Do metadata columns
meta_data_format <- srt_ctl_cpz_rcv@meta.data %>% rownames_to_column("cell") %>%
  left_join(meta_data %>% dplyr::select(cell, bio_rep, sample, type, annotation)) %>% 
  left_join(spagcn_map) %>%
  mutate(spagcn_clust = ifelse(is.na(spagcn_clust), "Unannotated", spagcn_clust)) %>%
  column_to_rownames("cell")

#Reinsert metadata
srt_ctl_cpz_rcv@meta.data <- meta_data_format

#Mito percentage calculate
srt_ctl_cpz_rcv <- PercentageFeatureSet(srt_ctl_cpz_rcv, "^mt-", col.name="percent_mito")

#Hemoglobin percentage calculate
srt_ctl_cpz_rcv <- PercentageFeatureSet(srt_ctl_cpz_rcv, "^Hb.*-", col.name="percent_hb")

### Do spot-level QC featureset plots 
p1 <- ggplot() +
  geom_histogram(data = srt_ctl_cpz_rcv@meta.data, aes(nFeature_Spatial), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Unique genes per spot") +
  geom_vline(xintercept = 100, linetype="dashed") +
  geom_vline(xintercept = 250, linetype="dashed") +
  geom_vline(xintercept = 500, linetype="dashed")

ggsave(plot = p1, filename=file.path(out_dir, "443_genes_per_spot.png"), width=5, height=5, units="in")

p1_tab <- srt_ctl_cpz_rcv@meta.data %>% 
  mutate(nFeature_Spatial_cut = case_when(nFeature_Spatial < 100 ~ "n<100",
                                          nFeature_Spatial >= 100 & nFeature_Spatial < 250 ~ "100<=n<250",
                                          nFeature_Spatial >= 250 & nFeature_Spatial < 500  ~ "250<=n<500",
                                          nFeature_Spatial >= 500 ~ "n>=500")
  ) %>% dplyr::count(bio_rep, nFeature_Spatial_cut) %>% pivot_wider(names_from = nFeature_Spatial_cut, values_from = n) %>%
  dplyr::select(bio_rep, `n<100`, `100<=n<250`, `250<=n<500`, `n>=500`)

fwrite(p1_tab, file = file.path(out_dir, "443_genes_per_spot_bio_rep_cut_tab.csv"))

p1_facet <- p1 + facet_wrap(~bio_rep, ncol = 3)

ggsave(plot = p1_facet, filename=file.path(out_dir, "443_genes_per_spot_bio_rep_facet.png"), width=5, height=5, units="in")



p1_overlay <- SpatialFeaturePlot(srt_ctl_cpz_rcv,features = "nFeature_Spatial",ncol=7) & 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits=c(min(srt_ctl_cpz_rcv@meta.data$nFeature_Spatial),max(srt_ctl_cpz_rcv@meta.data$nFeature_Spatial) ))

p2 <- ggplot() +
  geom_histogram(data = srt_ctl_cpz_rcv@meta.data, aes(nCount_Spatial), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Total counts per spots") 

ggsave(plot = p2, filename=file.path(out_dir, "443_counts_per_spot.png"), width=5, height=5, units="in")

p2_scatter <- ggplot() +
  geom_point(data = srt_ctl_cpz_rcv@meta.data, aes(x=nCount_Spatial, y=nFeature_Spatial), size=0.5) +
  geom_hline(yintercept = 500, linetype="dashed", color="red") +
  geom_vline(xintercept = 100, linetype="dashed", color="red")

ggsave(plot=p2_scatter, file = file.path(out_dir, "443_counts_per_spot_genes_per_spot_scatter.png"), width=5, height=5, units="in" )

p2_overlay <- SpatialFeaturePlot(srt_ctl_cpz_rcv,features = "nCount_Spatial",ncol=7) & 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits=c(min(srt_ctl_cpz_rcv@meta.data$nCount_Spatial),max(srt_ctl_cpz_rcv@meta.data$nCount_Spatial) ))

mito <- ggplot() +
  geom_histogram(data = srt_ctl_cpz_rcv@meta.data, aes(percent_mito), fill = "red", alpha = 0.7,  bins = 50) +
  ggtitle("Mitochondrial gene percentage") +
  geom_vline(xintercept=25, linetype="dashed") +
  geom_vline(xintercept=35, linetype="dashed") +
  geom_vline(xintercept=45, linetype="dashed")

ggsave(plot = mito, file.path(out_dir, "443_percent_mito_per_spot.png"), width=5, height=5, units="in" )

mito_facet <- mito + facet_wrap(~bio_rep, ncol = 3)

ggsave(plot = mito_facet, filename=file.path(out_dir, "443_percent_mito_per_spot_bio_rep_facet.png"), width=5, height=5, units="in")

mito_tab <- srt_ctl_cpz_rcv@meta.data %>% 
  mutate(percent_mito = case_when(percent_mito < 25 ~ "n<25",
                                  percent_mito >= 25 & percent_mito < 35 ~ "25<=n<35",
                                  percent_mito >= 35 & percent_mito < 45  ~ "35<=n<45",
                                  percent_mito >= 45 ~ "n>=45")
  ) %>% dplyr::count(bio_rep, percent_mito) %>% pivot_wider(names_from = percent_mito, values_from = n,values_fill = 0) %>%
  dplyr::select(bio_rep, `n<25`, `25<=n<35`, `35<=n<45`, `n>=45`)

fwrite(mito_tab, file = file.path(out_dir, "443_percent_mito_per_spot_bio_rep_cut_tab.csv"))

#Mito facet plot by annotation
mito_facet_annotation <- mito + facet_wrap(~annotation+bio_rep, ncol = 9)

ggsave(plot = mito_facet_annotation, filename=file.path(out_dir, "443_percent_mito_per_spot_bio_rep_annotation_facet.png"), width=1.5*9, height=1.5*6, units="in")

mito_tab_annotation <- srt_ctl_cpz_rcv@meta.data %>% 
  mutate(percent_mito = case_when(percent_mito < 25 ~ "n<25",
                                  percent_mito >= 25 & percent_mito < 35 ~ "25<=n<35",
                                  percent_mito >= 35 & percent_mito < 45  ~ "35<=n<45",
                                  percent_mito >= 45 ~ "n>=45")) %>% 
  dplyr::count(bio_rep, annotation, percent_mito) %>% 
  tidyr::complete(bio_rep, annotation, fill= list(n=0)) %>%
  pivot_wider(names_from = percent_mito, values_from = n,values_fill = 0) %>%
  dplyr::select(bio_rep, annotation, `n<25`, `25<=n<35`, `35<=n<45`, `n>=45`) %>% arrange(annotation, bio_rep)

fwrite(mito_tab_annotation, file = file.path(out_dir, "443_percent_mito_per_spot_bio_rep_annotation_cut_tab.csv"))

#Mito tab by spagcn cluster
mito_tab_spagcn_clust <- srt_ctl_cpz_rcv@meta.data %>% 
  mutate(percent_mito = case_when(percent_mito < 25 ~ "n<25",
                                  percent_mito >= 25 & percent_mito < 35 ~ "25<=n<35",
                                  percent_mito >= 35 & percent_mito < 45  ~ "35<=n<45",
                                  percent_mito >= 45 ~ "n>=45")) %>% 
  dplyr::count(bio_rep, spagcn_clust, percent_mito) %>% 
  tidyr::complete(bio_rep, spagcn_clust, fill= list(n=0)) %>%
  pivot_wider(names_from = percent_mito, values_from = n,values_fill = 0) %>%
  dplyr::select(bio_rep, spagcn_clust, `n<25`, `25<=n<35`, `35<=n<45`, `n>=45`) %>% arrange(spagcn_clust, bio_rep)

fwrite(mito_tab_spagcn_clust, file = file.path(out_dir, "443_percent_mito_per_spot_bio_rep_spagcn_clust_cut_tab.csv"))

p3_overlay <- SpatialFeaturePlot(srt_ctl_cpz_rcv,features = "percent_mito",ncol=7, combine=F) %>%
  lapply(function(x) x + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits=c(min(srt_ctl_cpz_rcv@meta.data$percent_mito),max(srt_ctl_cpz_rcv@meta.data$percent_mito) ))
  )

hb <- ggplot() +
  geom_histogram(data = srt_ctl_cpz_rcv@meta.data, aes(percent_hb), fill = "red", alpha = 0.7,  bins = 25) +
  ggtitle("Hemoglobin gene percentage") +
  geom_vline(xintercept=25, linetype="dashed") +
  geom_vline(xintercept=35, linetype="dashed") +
  geom_vline(xintercept=45, linetype="dashed")

#Note: there are 20 spots with percent_hb > 20, histogram is not useful due to binning
# srt_ctl_cpz_rcv@meta.data %>% filter(percent_hb > 20) %>% nrow
# [1] 20

hb_tab <- srt_ctl_cpz_rcv@meta.data %>% 
  mutate(percent_hb = case_when(percent_hb < 20 ~ "n<20",
                                percent_hb >= 20 & percent_hb < 30 ~ "20<=n<30",
                                percent_hb >= 30 & percent_hb < 40  ~ "30<=n<40",
                                percent_hb >= 40 ~ "n>=40")
  ) %>% dplyr::count(bio_rep, percent_hb) %>% pivot_wider(names_from = percent_hb, values_from = n,values_fill = 0) %>%
  dplyr::select(bio_rep, `n<20`, `20<=n<30`, `30<=n<40`, `n>=40`)

fwrite(hb_tab, file = file.path(out_dir, "443_percent_hb_per_spot_bio_rep_cut_tab.csv"))


### Make the gene-level QC plots
gene_attr <- data.frame(nUMI = Matrix::rowSums(srt_ctl_cpz_rcv@assays$Spatial@counts), 
                        nSpots = Matrix::rowSums(srt_ctl_cpz_rcv@assays$Spatial@counts > 0))

p3 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nUMI), fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10() +
  ggtitle("Total counts per gene (log10 scale)")

ggsave(plot=p3, file = file.path(out_dir, "443_counts_per_gene.png"), width=5, height=5, units="in" )
#Removed 7043 rows containing non-finite values (stat_bin)
# nrow(gene_attr %>% filter(nUMI == 0))
# [1] 7043

p4 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nSpots), fill = "red", alpha = 0.7,  bins = 50) +
  ggtitle("Total spots per gene")

ggsave(plot=p4, file = file.path(out_dir, "443_spots_per_gene.png"), width=5, height=5, units="in" )

p4_tab <- gene_attr %>% 
  mutate(nSpots_cut = case_when(nSpots < 6 ~ "n<6",
                                nSpots >= 6 & nSpots < 10 ~ "6<=n<10",
                                nSpots >= 10 & nSpots < 25  ~ "10<=n<25",
                                nSpots >= 25 ~ "n>=25")
  ) %>% dplyr::count(nSpots_cut) 

fwrite(p4_tab, file.path(out_dir, "443_spots_per_gene_cut_tab.csv"))

#Specific gene lookup
mt_genes <- grep("^mt-",rownames(srt_ctl_cpz_rcv@assays$Spatial@counts), value = T)
hb_genes <- grep("^Hb.*-",rownames(srt_ctl_cpz_rcv@assays$Spatial@counts), value = T)
ribo_genes <- grep("(^Rps|^Rpl)",rownames(srt_ctl_cpz_rcv@assays$Spatial@counts), value = T)
bc1_gene <- "Bc1"

#Nonzero genes
nz_genes <- gene_attr %>% filter(nUMI > 0)
                  
specific_gene_collapse <- map(lst(mt_genes,hb_genes,ribo_genes,bc1_gene),.f = ~paste(.,collapse=", "))
specific_gene_length <- map(lst(mt_genes,hb_genes,ribo_genes,bc1_gene), .f = length)

specific_gene_df <- setNames(stack(specific_gene_collapse),nm = c("gene names (collapsed)","category")) %>%
  left_join( setNames(stack(specific_gene_length),nm = c("gene number","category")) )

fwrite(specific_gene_df, file.path(out_dir, "443_specific_problematic_genes_tab.csv") )

#Get the top gene expressed plot
C = srt_ctl_cpz_rcv@assays$Spatial@counts
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]

png(filename = file.path(out_dir, "443_gene_top20_boxplot.png"), width = 15, height = 5, units = "in", res=300)
top20_boxplot <- boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, ylab = "gene % total count per cell",
        col = (scales::hue_pal())(20)[20:1])
mtext(top20_boxplot$names, at = 1:length(top20_boxplot$names), side = 1, line = 1, 
      col = case_when(top20_boxplot$names %in% mt_genes ~ "red",
                      top20_boxplot$names %in% hb_genes ~ "blue",
                      top20_boxplot$names %in% ribo_genes ~ "green",
                      top20_boxplot$names %in% bc1_gene ~ "purple",
                      TRUE ~ "black"
                      )
)
dev.off()


### Run the relevant QC steps

#Scenario 1 NBIS
  #Spot level 
  #>500 genes per spot
  #<25% mitochondrial gene reads percentage
  #<20% hemoglobin gene reads percentage
  
  #Gene level
  #Remove genes with all zeroes across all spots
  #Remove Bc1, Hb genes, mt genes

srt_scenario1 <- srt_ctl_cpz_rcv[, srt_ctl_cpz_rcv$nFeature_Spatial > 500 & 
                                   srt_ctl_cpz_rcv$percent_mito < 25 &
                                   srt_ctl_cpz_rcv$percent_hb < 20]

srt_scenario1 <- srt_scenario1[rownames(srt_scenario1) %in% rownames(nz_genes),]

gene_attr_anti <- gene_attr %>% rownames_to_column("gene_name") %>% filter(!gene_name %in% c(bc1_gene, hb_genes, mt_genes))

srt_scenario1 <- srt_scenario1[rownames(srt_scenario1) %in% gene_attr_anti$gene_name,]

#Save out the spot and gene csvs
srt_scenario1_spots <- data.frame("sampleid_spot" = colnames(srt_scenario1))
fwrite(srt_scenario1_spots, file.path(out_dir,"443_ST_QC_scenario1_spots.csv"))

srt_scenario1_genes <- data.frame("gene_name" = rownames(srt_scenario1))
fwrite(srt_scenario1_genes, file.path(out_dir,"443_ST_QC_scenario1_genes.csv"))

#Scenario 2 Edward Zhao
  #Spot level
  #>100 UMI per spot
  #<45% mitochondrial gene reads percentage
  
  #Gene level
  #Remove genes with all zeroes across all spots

srt_scenario2 <- srt_ctl_cpz_rcv[, srt_ctl_cpz_rcv$nCount_Spatial > 100 & 
                                   srt_ctl_cpz_rcv$percent_mito < 45]

#Save out the spot and gene csvs
srt_scenario2_spots <- data.frame("sampleid_spot" = colnames(srt_scenario2))
fwrite(srt_scenario2_spots, file.path(out_dir,"443_ST_QC_scenario2_spots.csv"))

srt_scenario2_genes <- data.frame("gene_name" = rownames(srt_scenario2))
fwrite(srt_scenario2_genes, file.path(out_dir,"443_ST_QC_scenario2_genes.csv"))

#Scenario 3 Hybrid
  #Spot level
  #>500 UMI per spot
  #<35% mitochondrial gene reads percentage
  
  #Gene level
  #Remove genes with all zeroes across all spots

srt_scenario3 <- srt_ctl_cpz_rcv[, srt_ctl_cpz_rcv$nFeature_Spatial > 500 & 
                                   srt_ctl_cpz_rcv$percent_mito < 35]

srt_scenario3 <- srt_scenario3[rownames(srt_scenario3) %in% rownames(nz_genes),]

#Save out the spot and gene csvs
srt_scenario3_spots <- data.frame("sampleid_spot" = colnames(srt_scenario3))
fwrite(srt_scenario3_spots, file.path(out_dir,"443_ST_QC_scenario3_spots.csv"))

srt_scenario3_genes <- data.frame("gene_name" = rownames(srt_scenario3))
fwrite(srt_scenario3_genes, file.path(out_dir,"443_ST_QC_scenario3_genes.csv"))

#Per contrast filters (would be applied at the contrast level)
  #Remove genes with <5 spots per gene (per DE contrast)

#Scenario 4 Sarbottam
#Spot level
#>500 UMI per spot
#<45% mitochondrial gene reads percentage

#Gene level
#Remove specific genes Bc1, Hb, Mt
#Remove genes with all zeroes across all spots

srt_scenario4 <- srt_ctl_cpz_rcv[, srt_ctl_cpz_rcv$nFeature_Spatial > 500 & 
                                   srt_ctl_cpz_rcv$percent_mito < 45]

srt_scenario4 <- srt_scenario4[rownames(srt_scenario4) %in% rownames(nz_genes),]

gene_attr_anti <- gene_attr %>% rownames_to_column("gene_name") %>% filter(!gene_name %in% c(bc1_gene, hb_genes, mt_genes))

srt_scenario4 <- srt_scenario4[rownames(srt_scenario4) %in% gene_attr_anti$gene_name,]

#Save out the spot and gene csvs
srt_scenario4_spots <- data.frame("sampleid_spot" = colnames(srt_scenario4))
fwrite(srt_scenario4_spots, file.path(out_dir,"443_ST_QC_scenario4_spots.csv"))

srt_scenario4_genes <- data.frame("gene_name" = rownames(srt_scenario4))
fwrite(srt_scenario4_genes, file.path(out_dir,"443_ST_QC_scenario4_genes.csv"))

#Scenario 5 Sarbottam
#Spot level
#>500 UMI per spot

#Gene level
#Remove genes with all zeroes across all spots

srt_scenario5 <- srt_ctl_cpz_rcv[, srt_ctl_cpz_rcv$nFeature_Spatial > 500]

srt_scenario5 <- srt_scenario5[rownames(srt_scenario5) %in% rownames(nz_genes),]

#Save out the spot and gene csvs
srt_scenario5_spots <- data.frame("sampleid_spot" = colnames(srt_scenario5))
fwrite(srt_scenario5_spots, file.path(out_dir,"443_ST_QC_scenario5_spots.csv"))

srt_scenario5_genes <- data.frame("gene_name" = rownames(srt_scenario5))
fwrite(srt_scenario5_genes, file.path(out_dir,"443_ST_QC_scenario5_genes.csv"))
