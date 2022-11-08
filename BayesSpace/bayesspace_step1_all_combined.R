library(scater)
library(BayesSpace)
library(viridis)
library(mclust)
library(harmony)
library(ggplot2)
setwd("/edgehpc/dept/compbio/projects/TST11523/mouse_cpz")
set.seed(100)

#Load data
sid = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sample.data = read.csv("data-raw/ID.map.csv") #sample metadata
sample.data$include='Y'
samples = unique(sample.data$include)
sample = samples[sid]
sample.ids = sample.data$SampleID[sample.data$include == sample]


sample.list = list() #reads data from Space Ranger output into R as a SingleCellExperiment
for(i in 1:length(sample.ids)){
  sample.id = sample.ids[i]
  sample.list[i] = readVisium(paste0("/edgehpc/dept/compbio/projects/TST11523/mouse_cpz/data-raw/", sample.id))
  
  sample.list[[i]]$SampleID = sample.id
}
names(sample.list) = sample.ids
sce.combined = do.call(cbind, sample.list) #merge into one object
rm(sample.list)

#QC
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce.combined)$gene_name)
sce.combined = sce.combined[!is_mito, ] #remove mitochondrial genes
sizeFactors = scuttle::librarySizeFactors(sce.combined)
sce.combined = sce.combined[,sizeFactors!=0] #remove spots with no expression


#Preprocess
if (sample == "Y"){
  sample.matrix = matrix(c(sample.ids), nrow = 7, ncol = 4) #offset spatial coordinates
} else if (sample == "CTL"){
  sample.matrix = matrix(c(sample.ids), nrow = 4, ncol = 3) #offset spatial coordinates
} else if (sample == "RCV"){
  sample.matrix = matrix(c(sample.ids), nrow = 3, ncol = 3) #offset spatial coordinates
}
sce.combined$row.orig = sce.combined$row
sce.combined$col.orig = sce.combined$col
for (i in 1:length(sample.ids)){
  coords = which(sample.matrix == sample.ids[i], arr.ind = T)
  sce.combined$row[sce.combined$SampleID == sample.ids[i]] = 80*coords[,"row"] + sce.combined$row.orig[sce.combined$SampleID == sample.ids[i]]
  sce.combined$col[sce.combined$SampleID == sample.ids[i]] = 130*coords[,"col"] + sce.combined$col.orig[sce.combined$SampleID == sample.ids[i]]
}
if(!dir.exists(paste0("out/bayes_combined/", sample))){
  dir.create(paste0("out/bayes_combined/", sample))
}
pdf(file = paste0("out/bayes_combined/", sample,"/capture.areas.", sample, ".pdf"), width = 10, height = 7)
clusterPlot(sce.combined, color = NA, "SampleID") +
  labs(fill = "Capture\narea")
dev.off()
sce.combined = spatialPreprocess(sce.combined, n.PCs = 50) #lognormalize, PCA
sce.combined = runUMAP(sce.combined, dimred = "PCA")
colnames(reducedDim(sce.combined, "UMAP")) = c("UMAP1", "UMAP2")

pdf(file = paste0("out/bayes_combined/", sample,"/umap.uncorrected.", sample, ".pdf"), width = 7, height = 7)
ggplot(data.frame(reducedDim(sce.combined, "UMAP")),
       aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$SampleID))) +
  geom_point(alpha = 0.1) +
  labs(color = "Sample") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw()
dev.off()

#Harmonization
sce.combined = RunHarmony(sce.combined, "SampleID", verbose = T)
sce.combined = runUMAP(sce.combined, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce.combined, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")
pdf(file = paste0("out/bayes_combined/", sample,"/umap.harmony.", sample, ".pdf"), width = 7, height = 7)
ggplot(data.frame(reducedDim(sce.combined, "UMAP.HARMONY")),
       aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$SampleID))) +
  geom_point(alpha = 0.1) +
  labs(color = "Sample") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw()
dev.off()

#BayesSpace tune clusters
mclust1 = Mclust(reducedDim(sce.combined, "HARMONY")[,1:20], G = 1:20, modelNames = "EEE")
pdf(file = paste0("out/bayes_combined/", sample,"/mclust.bic.", sample, ".pdf"), width = 10, height = 7)
plot(mclust1$BIC)
dev.off()

sce.combined = qTune(sce.combined, qs = seq(2,15), use.dimred = "HARMONY", nrep = 2000)
pdf(file = paste0("out/bayes_combined/", sample,"/bayesspace.qtune.", sample, ".pdf"), width = 10, height = 7)
qPlot(sce.combined)
dev.off()
saveRDS(sce.combined, paste0("data/bayes_combined/sce.combined.tuned.", sample, ".RDS"))
