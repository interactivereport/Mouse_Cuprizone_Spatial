library(tidyverse)
library(data.table)

#LR dataframe
lr_add_df <- data.frame(ligand=c("Il34","Apoe","Pecam1","Th","Th","Ddc","Ddc","Spp1"),
                        receptor=c("Csf1r","Trem2","Pdgfra","Drd1","Drd2","Drd1","Drd2","Cd44"),
                        species="Mouse")

out_dir <- file.path("/edgehpc/dept/compbio/projects/SingleCell/results/dev/spatial_tx/spatalk")

fwrite(lr_add_df, file.path(out_dir, "805_manual_lrs.csv"))

#Cell pairing df
man_sender_receiver <- data.frame(sender = c("microglia","microglia","microglia","astrocytes","mol","opc","mol","opc"),
                                  receiver = c("astrocytes","mol","opc","microglia","microglia","microglia","opc","mol"))

fwrite(man_sender_receiver, file.path(out_dir, "805_manual_sender_receiver.csv"))
