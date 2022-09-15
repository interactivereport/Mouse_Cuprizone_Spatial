library(tidyverse)
library(data.table)
library(glue)
library(Matrix)
library(optparse)
library(ggrepel)
library(clusterProfiler)

volcano_func <- function(result_df, FDR_threshold, FC_threshold){
  log2fc_threshold = log2(as.numeric(FC_threshold))
  res <- result_df %>% dplyr::filter(!is.na(FDR), FDR < 0.4) %>%
    dplyr::mutate(note = dplyr::case_when(
      FDR <= as.numeric(FDR_threshold) & log2FC >= log2(abs(as.numeric(FC_threshold))) ~ "Up-regulated",
      FDR <= as.numeric(FDR_threshold) & log2FC <= -log2(abs(as.numeric(FC_threshold))) ~ "Down-regulated",
      TRUE ~ "Not significant"
    )) %>%
    dplyr:::mutate(note = factor(note, levels = c("Up-regulated", "Down-regulated", "Not significant")))
  cols <- c("Up-regulated" = "#fc8d59", "Down-regulated" = "#91bfdb", "Not significant" = "#999999")
  volc <- ggplot(res, aes(y=-log10(FDR), x=log2FC, color=note)) +
    theme_bw() +
    theme(
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    geom_point(size=0.75) + 
    scale_colour_manual(values = cols) +
    geom_text_repel(data = res %>% filter(note %in% c("Up-regulated","Down-regulated")),
                    aes(label=ID), box.padding = 1, show.legend=FALSE) + 
    labs(x="log(FC)", y="-log10(Adj. p-value)",color="Threshold") +
    scale_x_continuous(expand=expansion(mult=c(0.02,0.02))) +
    scale_y_continuous(expand=expansion(mult=c(0.01,0.02)))
  return(volc)
}