library(tidyverse)
library(data.table)

out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/mouse_cuprizone_DE/010_geneset_celltypes_of_interest")

gsea_files <- list.files(out_dir, "gsea.csv", full.names=TRUE)

gsea_names <- gsub(".csv","",basename(gsea_files))

for(i in 1:length(gsea_files)){

gsea_df <- fread(gsea_files[[i]])

gsea_tab_format_top <- gsea_df %>% group_by(gs_name) %>%
  arrange(desc(NES)) %>% filter(p.adjust < 0.05) %>% slice(1:3)

gsea_tab_format_bottom <- gsea_df %>% group_by(gs_name) %>%
  arrange(NES) %>% filter(p.adjust < 0.05) %>% slice(1:3)

gsea_tab_format <- rbind(gsea_tab_format_top, gsea_tab_format_bottom) %>% distinct()

if(nrow(gsea_tab_format) > 0){
fwrite(gsea_tab_format,file.path(glue("/camhpc/home/mryals/edgehpc_biostat_results/{gsea_names[[i]]}_tab.csv")))
} else {message("skipped {gsea_names[[i]]}")}

}