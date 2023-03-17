
######################
# PCA plot of samples
# pseudobulk PCA
######################
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(Seurat)
library(ggbiplot,lib.loc = "/home/jluo/R_libs/")
library(openxlsx)
library(ggpubr)
library(svglite)
library(data.table)
library(sva)

# set directory
data_dir <- "/edgehpc/dept/compbio/users/ezhao/data-raw"
out_dir <- "/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx"

# metadata
samps <- read.xlsx(file.path(data_dir, "2020_07_22_Spatial_Transriptomics_Cuprizone_Mm_Pilot_Chip_Array_with_tissueAreas.xlsx"),
                   startRow = 5) %>% 
  dplyr::filter(Slide.ID != "") %>%
  mutate(samid=paste0(Slide.ID, Slide.Well)) %>% 
  column_to_rownames("samid")%>%
  mutate(label= gsub("0","", gsub("_","",Sample_ID)))

# # read individual rds files and create pseudobulk counts for each sample
# rds_dir <- file.path(data_dir, "rds")
# rds <- list.files(rds_dir, full.names = T)
# rds_list <- map(rds, readRDS)
# 
# # aggregate each row
# count <- map(rds_list, ~ rowSums(as.data.frame(.x@assays$Spatial@counts), na.rm = T))
# 
# # create count dataframe
# count_df <- as.data.frame(count) %>% 
#   `colnames<-`(str_replace(basename(rds),".RDS",""))
# fwrite(count_df, file.path(out_dir, "pseudobulk_raw_count.csv"))

# read it the once it was saved before, otherwise run the above code to generate the count_df object
count_df <- fread(file.path(out_dir, "pseudobulk_raw_count.csv"))

count_df <- count_df %>% 
  dplyr::select(rownames(samps))

################################################
# 0. No Batch_ID in the design model
################################################
# preliminary exploration
dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = samps,
                              design = ~ 1)
dds <- DESeq(dds)

vsd <- vst(dds, blind = T)
pcaData <- plotPCA(vsd, intgroup = c("TreatmentGroup"), returnData =TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData <- pcaData %>%
  left_join(samps %>% rownames_to_column("name") %>% dplyr::select(-c("TreatmentGroup","BioReplicate")),
            by=c("name"="name"))

g0 <- ggplot(pcaData, aes(x = PC1, y = PC2, shape=Batch_ID,
                    color = TreatmentGroup, group=TreatmentGroup)) +
  geom_point(size=0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_y_continuous(limits = c(-28, 20), expand = c(0,0))+
  scale_x_continuous(limits = c(-28, 20), expand = c(0,0))+
  # labs(title = "Batch 1&2 corrected")+
  # geom_text_repel(aes(label=label), show.legend=F)+
  stat_ellipse(lwd=line_size, color=1)+
  # geom_polygon(aes(color=TreatmentGroup))+
  coord_fixed()+
  theme_pubr()+
  theme(legend.position = "right",
        axis.text = element_text(size = font_size, family = ft_fmly),
        axis.title = element_text(size = font_size, family = ft_fmly),
        legend.text = element_text(size = font_size, family = ft_fmly),
        legend.title = element_text(size = font_size, family = ft_fmly),
        legend.key.size = unit(0.5, "line"),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size))

ggsave(file.path(out_dir, "pseudobulk_pca_wo_batch12_color_batch.png"),dpi = 600, width=3, height = 3)
ggsave(file.path(out_dir, "pseudobulk_pca_wo_batch12_color_batch.svg"), plot = g0, device = svg, width=3, height = 3)


######################################################################
# 1. Include Batch_ID in the design model to correct the batch effect
######################################################################

dds_batch <- DESeqDataSetFromMatrix(countData = count_df,
                       colData = samps,
                       design = ~ TreatmentGroup+Batch_ID)


dds_batch <- DESeq(dds_batch)
vsd <- vst(dds_batch, blind = F)
pcaData <- plotPCA(vsd, intgroup = c("TreatmentGroup", "Batch_ID"), returnData =TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData <- pcaData %>% 
  left_join(samps %>% rownames_to_column("name") %>% dplyr::select(-c("TreatmentGroup","BioReplicate","Batch_ID")), 
            by=c("name"="name")) 

font_size=7
line_size=0.3
ft_fmly="arial"
dds_batch <- DESeqDataSetFromMatrix(countData = count_df,
                                    colData = samps,
                                    design = ~ TreatmentGroup+Batch_ID)


dds_batch <- DESeq(dds_batch)
vsd <- vst(dds_batch, blind = F)
pcaData <- plotPCA(vsd, intgroup = c("TreatmentGroup", "Batch_ID"), returnData =TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData <- pcaData %>% 
  left_join(samps %>% rownames_to_column("name") %>% dplyr::select(-c("TreatmentGroup","BioReplicate","Batch_ID")), 
            by=c("name"="name")) 

font_size=7
line_size=0.3
ft_fmly="arial"

#Flip the color/shape aes
g_batch1_v2 <- ggplot(pcaData, aes(x = PC1, y = PC2, shape=Batch_ID,
                                   color = TreatmentGroup, group=TreatmentGroup)) +
  geom_point(size=0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_y_continuous(limits = c(-28, 20), expand = c(0,0))+
  scale_x_continuous(limits = c(-28, 20), expand = c(0,0))+
  # labs(title = "Batch 1&2 corrected")+
  # geom_text_repel(aes(label=label), show.legend=F)+
  stat_ellipse(lwd=line_size, color=1)+
  # geom_polygon(aes(color=TreatmentGroup))+
  coord_fixed()+
  theme_pubr()+
  theme(legend.position = "right",
        axis.text = element_text(size = font_size, family = ft_fmly),
        axis.title = element_text(size = font_size, family = ft_fmly),
        legend.text = element_text(size = font_size, family = ft_fmly),
        legend.title = element_text(size = font_size, family = ft_fmly),
        legend.key.size = unit(0.5, "line"),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size))

# ggsave(file.path(out_dir, "pseudobulk_pca_w_batch12_color_batch_v2.png"),dpi = 600, width=3, height = 3)
ggsave(file.path(out_dir, "pseudobulk_pca_w_batch12_color_batch_v2.svg"), plot = g_batch1_v2, device = svg, width=3, height = 3)


#######################################
# 2. batch correction using ComBat_seq
#######################################

corrected_samps <- samps %>% 
  mutate(BATCH=as.numeric(Batch_ID),
         GROUP=as.numeric(as.factor(TreatmentGroup)))
corrected_count_df <- ComBat_seq(count_df, group = corrected_samps$GROUP, batch = corrected_samps$BATCH) %>% 
  as.data.frame() %>% 
  `rownames<-`(rownames(count_df))



dds_batch_crt <- DESeqDataSetFromMatrix(countData = corrected_count_df,
                                    colData = samps,
                                    design = ~ TreatmentGroup)


dds_batch_crt <- DESeq(dds_batch_crt)
vsd <- vst(dds_batch_crt, blind = F)
pcaData <- plotPCA(vsd, intgroup = c("TreatmentGroup"), returnData =TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData <- pcaData %>% 
  left_join(samps %>% rownames_to_column("name") %>% dplyr::select(-c("TreatmentGroup","BioReplicate")), 
            by=c("name"="name")) 

font_size=7
line_size=0.3
ft_fmly="arial"

#Flip the color/shape aes
g_batch1_v2 <- ggplot(pcaData, aes(x = PC1, y = PC2, shape=Batch_ID,
                                   color = TreatmentGroup, group=TreatmentGroup)) +
  geom_point(size=0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_y_continuous(limits = c(-28, 20), expand = c(0,0))+
  scale_x_continuous(limits = c(-28, 22), expand = c(0,0))+
  # labs(title = "Batch 1&2 corrected")+
  # geom_text_repel(aes(label=label), show.legend=F)+
  stat_ellipse(lwd=line_size, color=1)+
  # geom_polygon(aes(color=TreatmentGroup))+
  coord_fixed()+
  theme_pubr()+
  theme(legend.position = "right",
        axis.text = element_text(size = font_size, family = ft_fmly),
        axis.title = element_text(size = font_size, family = ft_fmly),
        legend.text = element_text(size = font_size, family = ft_fmly),
        legend.title = element_text(size = font_size, family = ft_fmly),
        legend.key.size = unit(0.5, "line"),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size))
ggsave(file.path(out_dir, "pseudobulk_pca_w_batch12_corrected_color_batch_v1.png"),dpi = 600, width=3, height = 3)
ggsave(file.path(out_dir, "pseudobulk_pca_w_batch12_corrected_color_batch_v1.svg"), plot = g_batch1_v2, device = svg, width=3, height = 3)

