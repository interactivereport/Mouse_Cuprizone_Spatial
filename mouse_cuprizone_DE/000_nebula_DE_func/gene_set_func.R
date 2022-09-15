#Load packages
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(msigdbr)

#Do gene set analyses 
gene_set_func <- function(neb_df, ora_fdr, minGSSize){
  
  #Set a seed (seed also set inside gsea function)
  set.seed(2022)
  
  #Get gene sets (arbitrary)
  message(glue("Obtaining gene sets from msigdbr {packageVersion('msigdbr')}"))
  c2_pid = msigdbr(species = "mouse", category = "C2", subcategory = "CP:PID")
  c2_reactome = msigdbr(species = "mouse", category = "C2", subcategory = "CP:REACTOME")
  c2_wikipathways = msigdbr(species = "mouse", category = "C2", subcategory = "CP:WIKIPATHWAYS")
  c5_gobp = msigdbr(species = "mouse", category = "C5", subcategory = "GO:BP")
  
  gs_list <- list(c2_pid, c2_reactome, c2_wikipathways, c5_gobp)
  gs_names <- c("c2_pid","c2_reactome","c2_wikipathways","c5_gobp")
  
  #enricher
  neb_df_sig <- neb_df %>% dplyr::filter(FDR < ora_fdr) #No implied direction here
  enricher_res_list <- list()
  for(i in 1:length(gs_list)){
    msigdbr_df <- gs_list[[i]]
    msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
    if(nrow(neb_df_sig) > 0){
    enricher_obj <- enricher(gene = neb_df_sig$ID, 
                             TERM2GENE = msigdbr_t2g,
                             pAdjustMethod = "fdr",
                             minGSSize = minGSSize,
                             universe = neb_df$ID,
                             pvalueCutoff = 1)
    enricher_res_list[[i]] <- tryCatch({enricher_obj@result %>% arrange(p.adjust) }, error=function(e) NULL)
    } else{ enricher_res_list[[i]] <- NULL }
    
  }
  names(enricher_res_list) <- gs_names
  enricher_res_df <- rbindlist(enricher_res_list, idcol="gs_name")
  message("Finished ORA!")
  
  #fgsea
  neb_df_rank <- neb_df %>% select(ID, log2FC)
  neb_df_rank_vec <- neb_df_rank$log2FC
  names(neb_df_rank_vec) <- as.character(neb_df_rank$ID)
  neb_df_rank_vec <- sort(neb_df_rank_vec, decreasing=T)
  gsea_res_list <- list()
  for(i in 1:length(gs_list)){
    msigdbr_df <- gs_list[[i]]
    msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
    gsea_obj <- 
      GSEA(geneList=neb_df_rank_vec, 
                    TERM2GENE=msigdbr_t2g, 
                    minGSSize = minGSSize, eps = 0,
                    pvalueCutoff = 1, 
                    pAdjustMethod = "fdr", by = "fgsea", verbose=FALSE,seed = 1)
    gsea_res_list[[i]] <- tryCatch({ gsea_obj@result %>% arrange(p.adjust) }, error=function(e) NULL)
  }
  names(gsea_res_list) <- gs_names
  gsea_res_df <- rbindlist(gsea_res_list, idcol = "gs_name")
  message("Finished GSEA!")
  
  gene_set_res_list <- list("ora"=enricher_res_df,"gsea"=gsea_res_df)
  
  return(gene_set_res_list)
  
}