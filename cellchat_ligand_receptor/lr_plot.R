library(CellChat)
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(igraph)
library(ggpubr)
library(RColorBrewer)

# function to visualize the LR pairs (use cellchat_res only)
# to include not significant interactions, set thresh=1, otherwise 0.05
# to plot only interactions with p<0.01, set pval_cut=3
plot_lr <- function(cellchat=NULL, cellchat_res=NULL, emitter, receiver, thresh=1,
                    contra, cell_interest_role="source", 
                    LR_pair=NULL, cell_type_pair=NULL, pval_cut=NULL,regulated=TRUE,
                    source_target="target", cell_interest, direct_size=NULL){
  if (!is.null(cellchat)){
    idents <- unique(cellchat@idents$join)
    
    # extract data 
    # modify res data with thresh parameter changed to 1
    gg1 <- netVisual_bubble(cellchat, sources.use = idents,
                            targets.use = idents,  
                            comparison = c(1, 2), max.dataset = 2, 
                            title.name = "Increased in alternative group",
                            angle.x = 45, remove.isolate = T, thresh = thresh,
                            return.data=T)
    
    gg2 <- netVisual_bubble(cellchat, sources.use = idents,
                            targets.use = idents,  
                            comparison = c(1, 2), max.dataset = 1, 
                            title.name = "Decreased in alternative group", 
                            angle.x = 45, remove.isolate = T, thresh = thresh,
                            return.data=T)
    
    
    df1 <- gg1[["communication"]] %>% 
      data.frame() %>% 
      mutate(Regulated = "up")
    
    df2 <- gg2[["communication"]] %>% 
      data.frame() %>% 
      mutate(Regulated = "down") %>% 
      rbind(df1)

  
    
    #######################################################
    # create the rank plot
    df_rank <- df2 %>%
      # subset dataset based the cell type of interest and it's role (emitter or receiver)
      dplyr::filter(!!sym(cell_interest_role) == cell_interest) %>%
      group_by(dataset, !!sym(source_target)) %>%
      dplyr::summarise(N = n(),
                       N_sig=sum(pval > 1),
                       pct.sig=N_sig/N*100,
                       pb=sum(prob)) %>%
      dplyr::select(-N_sig) %>%
      mutate(dataset=case_when(dataset == "alt" ~ str_extract(contra, ".*?(?=vs)"),
                               dataset == "ref" ~ str_extract(contra, "(?<=vs).*"))) %>%
      pivot_wider(names_from = "dataset",values_from=c("pb","N", "pct.sig"), id_cols=source_target) %>%
      mutate(Delt_N=!!sym(paste0("pb_",str_extract(contra, ".*?(?=vs)")))-!!sym(paste0("pb_",str_extract(contra, "(?<=vs).*")))) %>%
      pivot_longer(cols = -c(all_of(source_target),Delt_N)) %>%
      separate(col = name, into = c("var", "dataset"), sep = "_") %>%
      pivot_wider(id_cols = c(all_of(source_target), Delt_N, dataset), names_from = "var", values_from = "value") %>%
      dplyr::filter(! is.na(.$N)) %>%
      mutate(pct_sig=ifelse(pct.sig <=20, 20,
                            ifelse(pct.sig <=40,  40,
                                   ifelse(pct.sig <=60, 60,
                                          ifelse(pct.sig <=80, 80, 100)))),
             N_cat=ifelse(N <=5, 5,
                          ifelse(N <=10,  10,
                                 ifelse(N <=15, 15,
                                        ifelse(N <=20, 20, 25)))),
             facet_lab=ifelse(is.na(Delt_N), "Missing in one", ifelse(Delt_N>0, str_extract(contra, ".*?(?=vs)"),str_extract(contra, "(?<=vs).*"))),
             unknow_cell=ifelse(str_detect(!!sym(source_target), "Unk"), 1, 0)) %>% 
      arrange(desc(Delt_N))
    
    if(cell_interest_role == "source"){role="Receiver"}else if(cell_interest_role=="target"){role="Emitter"}
    
    df_cont <- df_rank %>%
      dplyr::filter(facet_lab == str_extract(contra, "(?<=vs).*") | str_detect(facet_lab, "Missing")) %>%
      mutate(to_rm = ifelse(dataset == str_extract(contra, ".*?(?=vs)") & is.na(Delt_N), "yes","no")) %>%
      dplyr::filter(to_rm != "yes") %>%
      arrange(facet_lab, unknow_cell, Delt_N)
    
    value <- df_cont[source_target]
    levels =unique(df_cont[source_target]) %>%
      unlist()
    
    cont_p <- df_cont %>%
      mutate_at(vars(!!sym(source_target)), ~factor(.x, levels=rev(levels))) %>%
      ggplot(aes(x=!!sym(source_target), y=pb, color=dataset, alpha=pct_sig, size=N_cat))+
      geom_point()+
      labs(x= paste(role, "cell types"), y="Sum interaction probability",
           title = paste("Up-regulated in",str_extract(contra, "(?<=vs).*")))+
      coord_flip()+
      scale_color_discrete(name = "Groups")+
      scale_size_continuous(name = "No. of interactions", limits = c(5,100),
                            breaks=c(5,10,15,20,25),
                            labels=c("<5","5-10","10-15","15-20",">20"))+
      scale_alpha_continuous(name = "Pct of significant", limits=c(50,100), breaks=c(50,60,70, 80, 100))+
      theme(panel.background = element_blank(),
            legend.key = element_rect(fill = NA, color = NA),
            axis.line = element_line(colour = "black"),
            legend.position = "right",
            panel.grid.major.y = element_line(colour = "grey", linetype = 2),
            plot.background = element_rect(fill=NULL, linetype = NULL, color = NULL))
    
    df_alt <- df_rank %>%
      dplyr::filter(facet_lab == str_extract(contra, ".*?(?=vs)") | str_detect(facet_lab, "Missing")) %>%
      mutate(to_rm = ifelse(dataset == str_extract(contra, "(?<=vs).*") & is.na(Delt_N), "yes","no")) %>%
      dplyr::filter(to_rm != "yes") %>%
      arrange(facet_lab, unknow_cell, desc(Delt_N)) %>%
      mutate(dataset=="alt")
    
    value <- df_alt[source_target]
    levels =unique(df_alt[source_target]) %>%
      unlist()
    
    alt_p <- df_alt %>%
      mutate_at(vars(!!sym(source_target)), ~factor(.x, levels=rev(levels))) %>%
      ggplot(aes(x=!!sym(source_target), y=pb, color=dataset, alpha=pct_sig, size=N_cat))+
      geom_point()+
      labs(x= paste(role, "cell types"), y="Sum interaction probability",
           title = paste("Up-regulated in",str_extract(contra, ".*?(?=vs)")))+
      coord_flip()+
      scale_color_discrete(name = "Groups")+
      scale_size_continuous(name = "No. of interactions", limits = c(5,100),
                            breaks=c(5,10,15,20,25),
                            labels=c("<5","5-10","10-15","15-20",">20"))+
      scale_alpha_continuous(name = "Pct of significant", limits=c(50,100), breaks=c(50,60,70, 80, 100))+
      theme(panel.background = element_blank(),
            legend.key = element_rect(fill = NA, color = NA),
            axis.line = element_line(colour = "black"),
            legend.position = "right",
            panel.grid.major.y = element_line(colour = "grey", linetype = 2),
            plot.background = element_rect(fill=NULL, linetype = NULL, color = NULL))
    
    p_rank <- ggarrange(cont_p+rremove("xlab"), alt_p+rremove("xlab")+rremove("ylab"),
                        ncol = 2, common.legend = T, legend = "right")
    p_rank <- annotate_figure(p_rank, bottom = text_grob(label = "Sum interaction probability", hjust = 0.5, vjust = -0.25, rot = 0))+
      theme(plot.margin = unit(c(1,1,1,1), "cm"),
            plot.background = element_rect(fill=NULL, linetype = NULL, color = NULL))
    
    
    #################################################
    # create common interaction plot
    # prepare data for plotting common LR pairs
    
    
        
    if(is.null(pval_cut)){
      interact_union <- df2 %>% 
        arrange(source, target, interaction_name_2) %>% 
        select(-c(source.target, prob.original)) %>% 
        pivot_wider(names_from = "dataset", values_from = c("pval", "prob")) %>% 
        mutate(novel_ref=ifelse(is.na(pval_alt)&is.na(prob_alt),"novel_ref","NA"),
               novel_alt=ifelse(is.na(pval_ref)&is.na(prob_ref),"novel_alt","NA"),)     
    } else {
      interact_union <- df2 %>% 
        # dplyr::filter(pval==pval_cut) %>%
        mutate(Regulated=ifelse(pval<pval_cut, "NS", Regulated)) %>% 
        arrange(source, target, interaction_name_2) %>% 
        select(-c(source.target, prob.original)) %>% 
        pivot_wider(names_from = "dataset", values_from = c("pval", "prob")) %>% 
        mutate(novel_ref=ifelse(is.na(pval_alt)&is.na(prob_alt),1,0),
               novel_alt=ifelse(is.na(pval_ref)&is.na(prob_ref),1,0),)
    }
    
    
    # set missing interaction prob 0.01 and pval 1
    interact_union$prob_ref[is.na(interact_union$prob_ref)] <- 0.01
    interact_union$prob_alt[is.na(interact_union$prob_alt)] <- 0.01
    interact_union$pval_ref[is.na(interact_union$pval_ref)] <- 1
    interact_union$pval_alt[is.na(interact_union$pval_alt)] <- 1
    

    
    res <- interact_union %>% 
      pivot_longer(names_to = c("var","dataset"), names_sep = "_",
                   cols = c("prob_ref","prob_alt", "pval_ref","pval_alt", "novel_ref","novel_alt")) %>% 
      pivot_wider(names_from = "var", values_from = "value") %>% 
      mutate(source.target=paste0(source," -> ",target," (",dataset,")")) %>% 
      mutate(Regulated=ifelse(pval==1 & prob==0.01, "Not detected", Regulated)) %>% 
      mutate(interaction_name_flt=gsub("[[:space:]]", "", interaction_name_2)) %>% 
      dplyr::filter(!!sym(source_target) == cell_interest) %>% 
      mutate(up_down=factor(Regulated, levels = c("up","down","Not detected", "NS"))) %>%       
      mutate(plab=case_when(pval == 1 ~ "0.05",
                            pval == 2 ~ "0.01",
                            pval == 3 ~ "0.005"),
             prob_cat = ifelse(prob < 0.1, 0.1,
                               ifelse(prob < 0.5, 0.5,
                                      ifelse(prob < 1, 1,
                                             ifelse(prob < 1.5, 1.5, 2))))) %>% 
      mutate(Novel=ifelse(novel==1, 1, ""),
             Novel=factor(Novel, levels = c(1,0)))
    
 
    if (! is.null(cell_type_pair)){
      x_order <- cell_type_pair
    } else {
      x_order <- res %>% 
        setDT() %>% 
        .[, .N, by=c("group.names","dataset")] %>% 
        arrange(desc(N)) %>% 
        filter(dataset == "ref") %>% 
        pull("group.names")
      
    }
    
    ref_int <- res %>% 
      dplyr::filter(dataset=="ref") %>% 
      pull(interaction_name_2)
    alt_int <- res %>% 
      dplyr::filter(dataset=="alt") %>% 
      pull(interaction_name_2)
    common <- intersect(ref_int, alt_int)
    
    df_common <- res %>% 
      dplyr::filter(interaction_name_2 %in% common)
    
    ref_df <- df_common %>% 
      dplyr::filter(dataset == "ref") %>%
      dplyr::filter(group.names %in% x_order) %>% 
      mutate(source.target = gsub("ref", str_extract(contra, "(?<=vs).*"), source.target),
             group.names=factor(group.names, levels = x_order))
    
    # unique(df_common$group.names)
    alt_df <- df_common %>% 
      dplyr::filter(dataset == "alt") %>% 
      dplyr::filter(group.names %in% x_order) %>%       
      mutate(source.target = gsub("alt", str_extract(contra, ".*?(?=vs)"), source.target),
             group.names=factor(group.names, levels = x_order)) 
  }else{print("Please provide the cellchat object!")}
  
  # select the top 15 sum of prob LR pairs OR user specified LR for plotting
  if (! is.null(LR_pair)){
    ord1 <- LR_pair[LR_pair %in% ref_df$interaction_name_flt]  
  } else {
    ord1 <- ref_df %>% 
      setDT() %>% 
      .[,lapply(.SD, sum, na.rm=T), .SD=c("prob"), by=interaction_name_flt] %>% 
      arrange(desc(prob)) %>% 
      top_n(15) %>% 
      pull(interaction_name_flt)
  }
  # order the LR pairs based on the number of cell-type with this interaction
  y_order_ref <- ref_df %>% 
    dplyr::filter(interaction_name_flt %in% ord1) %>% 
    setDT() %>% 
    .[, .N, by=c("interaction_name_2")] %>% 
    arrange(N) %>% 
    pull("interaction_name_2")
  
  if (regulated){
    ref_df <- ref_df %>%
      dplyr::filter(interaction_name_flt %in% ord1) %>%
      mutate(interaction_name_2=factor(interaction_name_2, levels = y_order_ref))
    
    ref_g <- ref_df %>%
      ggplot(aes(x=group.names, y=interaction_name_2, color=plab, size=prob_cat))+
      geom_point()+
      geom_point(aes(shape=as.factor(up_down)), fill="black", size=direct_size)+
      labs(x="Emitter-receiver pair", y="Ligand-receptor pair", title = str_extract(contra, "(?<=vs).*"))+
      # scale_color_discrete(name = "Significance level", breaks=c("0.005", "0.01", "0.05"),
      #                      labels=c("<0.01", "0.05-0.01", ">0.05"))+
      scale_color_manual(name = "Significance level",
                         values=c("0.005"="#F8766D", "0.01"="#00BA38", "0.05"="#619CFF"),
                         labels=c("<0.01", "0.05-0.01", ">0.05"),drop=F)+
      scale_size_continuous(name = "Interaction probability", limits = c(0.1,2),
                            breaks=c(0.1, 0.5, 1, 1.5, 2),
                            labels = c("<0.1","0.1-0.5","0.5-1.0","1.0-1.5",">1.5"))+
      scale_shape_manual(name="up_down", values = c(24,25,4,8),
                         labels=c("up","down","Not detected","NS"))+
      theme(panel.background = element_blank(),
            legend.key = element_rect(fill = NA, color = NA),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            panel.grid.major.y = element_line(colour = "grey", linetype = 2))
  }else{
    ref_df1 <- ref_df %>% 
      dplyr::filter(interaction_name_flt %in% ord1) %>% 
      mutate(interaction_name_2=factor(interaction_name_2, levels = y_order_ref)) %>% 
      select(source, target, group.names, interaction_name_2, up_down, prob_cat, Regulated, interaction_name_flt, Novel)
    
    # for LR not present
    lrs1 <- LR_pair[!LR_pair %in% ref_df$interaction_name_flt]
    
    ref_df2 <- data.frame(group.names=rep(x_order, length(lrs1)),
                          interaction_name_flt=rep(lrs1, each=length(x_order))) %>% 
      mutate(source=str_extract(group.names, ".*?(?= )"),
             target=str_extract(group.names, "(?<= -> ).*"),
             interaction_name_2=gsub("-", " - ", interaction_name_flt),
             up_down="Not detected",
             prob_cat=0.5,
             Regulated="Not detected",
             Novel=NA) %>%
      rbind(ref_df1) %>% 
      mutate(dataset=str_extract(contra, "(?<=vs).*"))
    
    pval_cut_value <- ifelse(pval_cut==3, 0.01, 0.05)
    ref_g <- ref_df2 %>% 
      ggplot(aes(x=group.names, y=interaction_name_flt, color=Regulated, size=prob_cat))+
      geom_point()+
      geom_point(aes(shape=Novel,size=prob_cat*2))+
      labs(x="Emitter-receiver pair", y="Ligand-receptor pair", title = str_extract(contra, "(?<=vs).*"))+
      scale_color_manual(name = "up_down",
                         values=c("up"="#fc8d59", "down"="#91bfdb", "Not detected"="grey", "NS"="black"),
                         labels=c("up","down","Not detected",paste0("NS (p<",pval_cut_value,")")),
                         drop=F)+
      scale_size_continuous(name = "Interaction probability", limits = c(0.1,2),
                            breaks=c(0.1, 0.5, 1, 1.5, 2),
                            labels = c("<0.1","0.1-0.5","0.5-1.0","1.0-1.5",">1.5"))+
      scale_shape_manual(name="Novel", values = c(8,NULL,NULL),
                         labels=c("novel","",""))+
      theme(panel.background = element_blank(),
            legend.key = element_rect(fill = NA, color = NA),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            panel.grid.major.y = element_line(colour = "grey", linetype = 2),
            plot.background = element_rect(fill=NULL, linetype = NULL, color = NULL))    
  }
  
  
  if (! is.null(LR_pair)){
    ord2 <- LR_pair[LR_pair %in% alt_df$interaction_name_flt]
  } else {
    ord2 <- alt_df %>% 
      setDT() %>% 
      .[,lapply(.SD, sum, na.rm=T), .SD=c("prob"), by=interaction_name_flt] %>% 
      arrange(desc(prob)) %>% 
      top_n(15) %>% 
      pull(interaction_name_flt)
  }
  
  # order the LR pairs based on the number of cell-type with this interaction  
  y_order_alt <- alt_df %>% 
    dplyr::filter(interaction_name_flt %in% ord2) %>% 
    setDT() %>% 
    .[, .N, by=c("interaction_name_2")] %>% 
    arrange(N) %>% 
    pull("interaction_name_2")
  
  if (regulated){
    alt_df <- alt_df %>% 
      dplyr::filter(interaction_name_flt %in% ord2) %>% 
      mutate(interaction_name_2=factor(interaction_name_2, levels = y_order_alt)) 
    alt_g <- alt_df %>%     
      ggplot(aes(x=group.names, y=interaction_name_2, color=plab, size=prob_cat))+
      geom_point()+
      geom_point(aes(shape=as.factor(up_down)), fill="black", size=direct_size)+
      labs(x="Emitter-receiver pair", y="Ligand-receptor pair", title = str_extract(contra, ".*?(?=vs)"))+
      scale_color_manual(name = "Significance level", 
                         values=c("0.005"="#F8766D", "0.01"="#00BA38", "0.05"="#619CFF"),
                         labels=c("<0.01", "0.05-0.01", ">0.05"),drop=F)+
      scale_size_continuous(name = "Interaction probability", limits = c(0.1,2), 
                            breaks=c(0.1, 0.5, 1, 1.5, 2), 
                            labels = c("<0.1","0.1-0.5","0.5-1.0","1.0-1.5",">1.5"))+
      scale_shape_manual(name="Regulated", values = c(24,25,4,8), 
                         labels=c("up","down","Not detected","NS"))+
      theme(panel.background = element_blank(),
            legend.key = element_rect(fill = NA, color = NA),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            panel.grid.major.y = element_line(colour = "grey", linetype = 2))    
  }else{
    alt_df1 <- alt_df %>% 
      dplyr::filter(interaction_name_flt %in% ord2) %>% 
      mutate(interaction_name_2=factor(interaction_name_2, levels = y_order_alt)) %>% 
      select(source, target, group.names, interaction_name_2, up_down, prob_cat, Regulated, interaction_name_flt, Novel)
    
    # for LR not present
    lrs2 <- LR_pair[!LR_pair %in% alt_df$interaction_name_flt]
    
    alt_df2 <- data.frame(group.names=rep(x_order, length(lrs2)),
                          interaction_name_flt=rep(lrs2, each=length(x_order))) %>% 
      mutate(source=str_extract(group.names, ".*?(?= )"),
             target=str_extract(group.names, "(?<= -> ).*"),
             interaction_name_2=gsub("-", " - ", interaction_name_flt),
             up_down="Not detected",
             prob_cat=0.5,
             Regulated="Not detected",
             Novel=NA) %>%
      rbind(alt_df1) %>% 
      mutate(dataset=str_extract(contra, ".*?(?=vs)"))

    
    pval_cut_value <- ifelse(pval_cut==3, 0.01, 0.05)
    alt_g <-
      alt_df2 %>%     
      ggplot(aes(x=group.names, y=interaction_name_flt, color=Regulated, size=prob_cat))+
      geom_point()+
      geom_point(aes(shape=Novel,size=prob_cat*2))+
      labs(x="Emitter-receiver pair", y="Ligand-receptor pair", title = str_extract(contra, ".*?(?=vs)"))+
      scale_color_manual(name = "up_down",
                         values=c("up"="#fc8d59", "down"="#91bfdb", "Not detected"="grey", "NS"="black"),
                         labels=c("up","down","Not detected",paste0("NS (p<",pval_cut_value,")")),
                         drop=F)+
      scale_size_continuous(name = "Interaction probability", limits = c(0.1,2), 
                            breaks=c(0.1, 0.5, 1, 1.5, 2), 
                            labels = c("<0.1","0.1-0.5","0.5-1.0","1.0-1.5",">1.5"))+
      scale_shape_manual(name="Novel", values = c(8,NULL,NULL),
                         labels=c("novel","",""))+     
      theme(panel.background = element_blank(),
            legend.key = element_rect(fill = NA, color = NA),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            panel.grid.major.y = element_line(colour = "grey", linetype = 2),
            plot.background = element_rect(fill=NULL, linetype = NULL, color = NULL))    
  }
  
  
  g <- ggarrange(ref_g+rremove("xlab")+rremove("ylab"), alt_g+rremove("ylab"), 
                 ncol = 1, common.legend = T, legend = "right")
  g <- annotate_figure(g, left = text_grob(label = "Ligand-receptor pair", hjust = 0.5, vjust = 0.5, rot = 90))+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.background = element_rect(fill=NULL, linetype = NULL, color = NULL))
  
  # rank_plot = p_rank, rank_data = df_rank,
  res <- list(rank_plot = p_rank, rank_data = df_rank, 
              common_plot = g, common_data=rbind(ref_df2, alt_df2))
  return(res)
  
}

