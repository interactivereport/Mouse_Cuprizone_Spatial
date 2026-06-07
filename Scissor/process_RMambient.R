# yarrow06 (1TB memory)
# use conda (/home/zouyang/anaconda3/): conda activate /home/zouyang/anaconda3/envs/Scissor
# export OPENBLAS_NUM_THREADS=1
# nohup Rscript process_RMambient.R > result/region_RMambient/log &
# vignetes: https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html
#conda activate Scissor
#export OPENBLAS_NUM_THREADS=1
#nohup Rscript process.R > result/region_RMambient/log &

## loading ------------
rm(list=ls())
gc()
loadPKG <- function(){
  require(rhdf5)
  require(Matrix)
  require(Seurat)
  require(Scissor)
  source("Scissor.R")
}
suppressMessages(suppressWarnings(loadPKG()))
getOBS <- function(strH5ad){
  message("\tobtainning obs ...")
  obs <- h5read(strH5ad,"obs")
  sel <- names(obs)[sapply(obs,function(x)return(is.null(names(x))))&!grepl("^_|index$",names(obs))]
  meta <- do.call(cbind.data.frame, obs[sel])
  #meta <- do.call(cbind.data.frame, obs[grep("^_",names(obs),invert=T)])
  #dimnames(meta) <- list(obs[["_index"]],grep("^_",names(obs),invert=T,value=T))
  rownames(meta) <- obs[[grep("index$",names(obs))]]
  for(one in names(obs[["__categories"]])){
    if(min(meta[,one])<0){
      ann <- meta[,one]+1
      ann[ann<1] <- max(ann)+1
      annLable <- obs[["__categories"]][[one]]
      annLable <- c(annLable,"NAN")
      meta[,one] <- annLable[ann]
    }else{
      meta[,one] <- obs[["__categories"]][[one]][1+meta[,one]]
    }
  }
  # for anndata v 0.8
  for(one in names(obs)[sapply(obs,function(x)return(!is.null(names(x))))&!grepl("^_|index$",names(obs))]){
    if(sum(c("categories","codes")%in%names(obs[[one]]))!=2) next
    meta[[one]] <- obs[[one]]$categories[1+obs[[one]]$codes]
  }
  return(meta)
}
as_matrix <- function(SM,v=0L){
  message("\t transform to the dense matrix")
  rN <- nrow(SM)
  cN <- ncol(SM)
  tmp <- matrix(v,nrow=rN,ncol=cN,dimnames = dimnames(SM))
  if(is.integer(v)) SM@x <- ifelse(SM@x!=0,1L,0L)
  
  cI <- findInterval(seq(SM@x)-1,SM@p[-1])+1
  for(i in seq_along(SM@x))
    tmp[1+SM@i[i],cI[i]] <- SM@x[i]
  return(tmp)
}
getobsm <- function(strH5ad,key){
  AllK <- h5ls(strH5ad,recursive=2)
  k <- AllK[grepl("obsm",AllK[,1]),2]
  if(!key%in%k) return(NULL)
  X <- h5read(strH5ad,paste0("obsm/",key))
  colnames(X) <- getID(strH5ad,AllK,"/obs")    #h5read(strH5ad,"/obs/_index")
  return(t(X))
}
getobsmKey <- function(strH5ad){
  k <- h5ls(strH5ad,recursive=2)
  k <- k[grepl("obsm",k[,1]),2]
  return(k)
}
getID <- function(strH5ad,keys,grp){
  if("_index" %in% keys$name[grepl(grp,keys$group)]){
    return(h5read(strH5ad,paste0(grp,"/_index")))
  }else if("index" %in% keys$name[grepl(grp,keys$group)]){
    return(h5read(strH5ad,paste0(grp,"/index")))
  }else{
    stop(paste("unknown adata format: Neither index or _index exists in group",grp))
  }
}
#
strH5ad <- "/camhpc/ngs/projects/TST11621/Scissor/data/cuprizonemouse_annotated_Ambient_RNA_Doublet_removed.h5ad"
strBulk <- "/camhpc/ngs/projects/TST11391/dnanexus/20220525205955_Wei.Li.EXTERNAL./"
strOut <- "result/region_RMambient/"

## Bulk (TPM) and phenotype ---------
message("Bulk: get TPM")
TPM <- data.table::fread(paste0(strBulk,"genes.tpm.tsv"))
TPM <- as.matrix(data.frame(row.names=unlist(TPM[,1]),TPM[,-1],check.names = F))
gInfo <- data.table::fread("/camhpc/ngs/genomes/DNAnexus_references/rnaseq/mouse/Mouse.GRCm38.vM25.l1_5.ERCC/Mouse.GRCm38.vM25.l1_5.ERCC.transcript.gene_info.csv")
gInfo <- data.frame(row.names=unlist(gInfo$geneID),gInfo)
TPM <- TPM[rownames(gInfo),]

# merge transcripts of the same gene name
for(one in unique(gInfo$gene_name[duplicated(gInfo$gene_name)])){
  if(nchar(one)<1) next
  ix <- which(gInfo$gene_name%in%one)
  TPM[ix[1],] <- apply(TPM[ix,],2,sum)
}
TPM <- TPM[!duplicated(gInfo$gene_name)&nchar(gInfo$gene_name)>1,]
rownames(TPM) <- unlist(gInfo[rownames(TPM),"gene_name"])

## sc: create seurat object from h5ad ------
message("SC: get cell annotation")
meta <- getOBS(strH5ad)
cID <- rownames(meta)
strSub <- paste0(strOut,"selectedCell.rds")
if(file.exists(strSub)){
  selC <- readRDS(strSub)
}else{
  selC <- cID#sample(cID,ceiling(length(cID)/2))#
  #saveRDS(selC,strSub)
}
meta <- meta[selC,]
for(oneK in getobsmKey(strH5ad)){
  X <- getobsm(strH5ad,oneK)
  colnames(X) <- paste(gsub("^X_","",oneK),1:ncol(X),sep="_")
  meta <- merge(meta,X,by="row.names",all.x=TRUE)
  rownames(meta) <- meta[,1]
  meta <- meta[,-1]
}
saveRDS(meta,file=paste0(strOut,"meta.rds"))
# get gene information
message("SC: get gene information")
gInfo <- h5read(strH5ad,"raw/var")
gID <- gInfo[[grep("index$",names(gInfo),value=T)[1]]]
selG <- intersect(rownames(TPM),gID)

## get counts
message("SC: get counts")
X <- h5read(strH5ad,"raw/X")
X <- sparseMatrix(i=X$indices+1,p=X$indptr,x=as.numeric(X$data),
                  dims=c(length(gID),length(cID)),
                  dimnames=list(gID,cID))
X <- X[gID%in%selG,cID%in%selC]
scaleF <- as.numeric(quantile(colSums(X),0.95))

## create seurat object
message("SC: create seurat object")
#sc_dataset <- Seurat_preprocessing(X, verbose = F)
#suppressWarnings({sc_dataset <- AddMetaData(sc_dataset,meta)})
suppressWarnings({
  sc_dataset <- CreateSeuratObject(X,
                                   project="TST11621_Scissor",
                                   min.cells=0,min.features=0,
                                   meta.data=meta)
})
rm(X)
gc()
sc_dataset <- NormalizeData(sc_dataset,
                            normalization.method = "LogNormalize",
                            scale.factor = scaleF)

# get the embedding layout
message("SC: get layout embedding")
suppressWarnings({
  for(one in c("pca","umap")){#,"umap_celltype4","umap_celltype6"
    embedd <- matrix(t(h5read(strH5ad,paste0("obsm/X_",one))),nrow=length(cID),dimnames=list(cID,NULL))[selC,]
    sc_dataset[[one]] <- CreateDimReducObject(embeddings=embedd,key=one,assay="RNA")
  }
})

# create kNN and sNN from PC of h5ad
message("SC: Seurat FindNeighbors according to PCA from h5ad")
sc_dataset <- FindNeighbors(sc_dataset,dims=1:50)
strCluster <- paste0(strOut,"cluster.rds")
if(!file.exists(strCluster)){
  res <- 0.001
  sc_id <- as.factor(paste(sc_dataset@meta.data$Region,sc_dataset@active.ident,sep="_"))
  #browser()
  while(max(table(sc_id))>160000){
    message("resolution: ",res)
    sc_dataset <- FindClusters(sc_dataset,resolution=res)
    sc_id <- as.factor(paste(sc_dataset@meta.data$Region,sc_dataset@active.ident,sep="_"))
    res <- res+0.0001
  }
  names(sc_id) <- names(sc_dataset@active.ident)
  saveRDS(sc_id,file=strCluster)
}else{
  sc_id <- readRDS(strCluster)
}

# saving memory remove unnecessary entries and obtain the dense matrix of sNN
sc_data <- as_matrix(sc_dataset@assays$RNA@data,0)
sc_snn <- sc_dataset@graphs$RNA_snn
rm(sc_dataset)
diag(sc_snn) <- 0

## scissor: TPM, D ---------
message("Scissor:")
meta <- gsub("TPM\\|","",sapply(strsplit(colnames(TPM),"\\-"),head,1))

for(oneP in list(c("cpz4W","control"),c("Veh","cpz4W"))){
  message("\t",paste(oneP,collapse=".vs."))
  Bulk <- TPM[,meta%in%oneP]
  phenotype <- ifelse(gsub("TPM\\|","",sapply(strsplit(colnames(Bulk),"\\-"),head,1))==oneP[1],1,0)
  
  for(i in levels(sc_id)){
    selC <- names(sc_id)[sc_id==i]
    message("\n\n\t*** Cluster ",i,": ",length(selC)," cells")
    if(length(selC)==0) next
    sc_sub_snn <- as_matrix(sc_snn[cID%in%selC,cID%in%selC])
    sciCell <- Scissor_custom(Bulk,
                              sc_data[,cID%in%selC],sc_sub_snn,
                              phenotype, tag=rev(oneP),
                              alpha = 0.05, family = "binomial",
                              Save_file = paste0(strOut,'Scissor_',paste(oneP,collapse=".vs."),"_cluster",i,'.RData'))
    saveRDS(sciCell,paste0(strOut,'Scissor_',paste(oneP,collapse=".vs."),"_cluster",i,'.rds'))
  }
}

## -----
