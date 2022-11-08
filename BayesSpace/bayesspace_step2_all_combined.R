library(BayesSpace)

setwd("/edgehpc/dept/compbio/projects/TST11523/mouse_cpz")
set.seed(100)

#Load data
sid = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sample.data = read.csv("data-raw/ID.map.csv") #sample metadata

sample.data$include='Y'
types = unique(sample.data$include)

qs = c(14,15,16,17,18,19,20)
array = expand.grid(q = qs, Type = types, stringsAsFactors = FALSE) # sample = samples, 
q = array$q[sid]
#sample = array$sample[sid]
type = array$Type[sid]
sce.combined = readRDS(paste0("data/bayes_combined/sce.combined.tuned.", type, ".RDS"))

sce.combined = spatialCluster(sce.combined, q = q,
                              use.dimred = "HARMONY", nrep = 25000)
outRDS = paste0("data/bayes_combined/sce.bayesspace.", type, ".", q, ".RDS")
saveRDS(sce.combined, file = outRDS)

pdf(file=paste0("out/bayes_combined/", sample, "/sce.bayesspace.", type, ".", q, ".pdf"),
    width = 18, height = 10)
clusterPlot(sce.combined, color = NA)
dev.off()
