#Assign seed inside generate newmeta parallel
.generate_newmeta_doParallel_v2 <- function (st_meta, st_dist, min_percent, rand_seed) {
  st_meta <- st_meta[st_meta$label != "less nFeatures", ]
  cellname <- colnames(st_meta)[-c(1:7)]
  newmeta <- foreach::foreach(i = 1:nrow(st_meta), .combine = rbind, 
                              .packages = "Matrix", .export = c(".det_neighbor", ".get_weight1", 
                                                                ".get_weight2")) %dopar% {
                                                                  set.seed(rand_seed)
                                                                  newmeta_spot <- NULL
                                                                  newmeta_ratio <- NULL
                                                                  newmeta_cell <- NULL
                                                                  newmeta_x <- NULL
                                                                  newmeta_y <- NULL
                                                                  spot_name <- st_meta$spot[i]
                                                                  spot_x <- st_meta$x[i]
                                                                  spot_y <- st_meta$y[i]
                                                                  spot_cellnum <- st_meta$cell_num[i]
                                                                  spot_percent <- as.numeric(st_meta[i, -c(1:7)])
                                                                  spot_percent <- spot_percent * spot_cellnum
                                                                  spot_percent_floor <- floor(spot_percent)
                                                                  spot_percent_dec <- spot_percent - spot_percent_floor
                                                                  spot_celltype <- which(spot_percent_dec > min_percent)
                                                                  k <- 0
                                                                  if (length(spot_celltype) > 0) {
                                                                    spot_percent <- spot_percent_dec[spot_celltype]
                                                                    spot_celltype <- cellname[spot_celltype]
                                                                    for (j in 1:length(spot_celltype)) {
                                                                      spot_percent1 <- spot_percent[j]
                                                                      newmeta_ratio <- c(newmeta_ratio, spot_percent1)
                                                                      newmeta_cell <- c(newmeta_cell, spot_celltype[j])
                                                                      k <- k + 1
                                                                    }
                                                                  }
                                                                  spot_celltype <- which(spot_percent_floor > 0)
                                                                  if (length(spot_celltype) > 0) {
                                                                    spot_percent <- spot_percent_floor[spot_celltype]
                                                                    spot_celltype <- cellname[spot_celltype]
                                                                    spot_cell <- NULL
                                                                    for (j in 1:length(spot_celltype)) {
                                                                      spot_cell <- c(spot_cell, rep(spot_celltype[j], 
                                                                                                    spot_percent[j]))
                                                                    }
                                                                    for (j in 1:length(spot_cell)) {
                                                                      spot_percent1 <- 1
                                                                      newmeta_ratio <- c(newmeta_ratio, spot_percent1)
                                                                      newmeta_cell <- c(newmeta_cell, spot_cell[j])
                                                                      k <- k + 1
                                                                    }
                                                                  }
                                                                  if (k > 0) {
                                                                    newmeta_spot <- c(newmeta_spot, rep(spot_name, k))
                                                                    if (k == 1) {
                                                                      newmeta_x <- c(newmeta_x, spot_x)
                                                                      newmeta_y <- c(newmeta_y, spot_y)
                                                                    }
                                                                    else {
                                                                      n_neighbor <- 4
                                                                      st_dist1 <- st_dist[, spot_name]
                                                                      st_dist1 <- st_dist1[st_dist1 > 0]
                                                                      st_dist1 <- st_dist1[order(st_dist1)]
                                                                      st_dist1 <- st_dist1[1:n_neighbor]
                                                                      st_meta_neighbor <- st_meta[st_meta$spot %in% 
                                                                                                    names(st_dist1), ]
                                                                      st_meta_neighbor <- .det_neighbor(st_meta_neighbor, 
                                                                                                        spot_x, spot_y, st_dist1)
                                                                      if (nrow(st_meta_neighbor) == 0) {
                                                                        for (j in 1:k) {
                                                                          st_angle_new <- sample(x = c(1:360), size = 1)
                                                                          st_dist_new <- sample(x = c(0:st_dist1), 
                                                                                                size = 1)
                                                                          newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * 
                                                                                                                     pi/180)
                                                                          newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * 
                                                                                                                     pi/180)
                                                                          newmeta_x <- c(newmeta_x, newmeta_x1)
                                                                          newmeta_y <- c(newmeta_y, newmeta_y1)
                                                                        }
                                                                      }
                                                                      else {
                                                                        for (j in 1:k) {
                                                                          sc_name <- newmeta_cell[j]
                                                                          sc_w1 <- .get_weight1(st_meta_neighbor, sc_name)
                                                                          set.seed(j)
                                                                          st_angle_new <- sample(x = c(1:360), size = 1, 
                                                                                                 prob = as.numeric(sc_w1))
                                                                          spot_ratio <- st_meta[st_meta$spot == spot_name, 
                                                                                                sc_name]
                                                                          st_dist_new <- .get_weight2(st_meta_neighbor, 
                                                                                                      sc_name, st_dist1, st_angle_new, spot_ratio)
                                                                          newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * 
                                                                                                                     pi/180)
                                                                          newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * 
                                                                                                                     pi/180)
                                                                          newmeta_x <- c(newmeta_x, newmeta_x1)
                                                                          newmeta_y <- c(newmeta_y, newmeta_y1)
                                                                        }
                                                                      }
                                                                    }
                                                                    data.frame(spot = newmeta_spot, cell_ratio = newmeta_ratio, 
                                                                               celltype = newmeta_cell, x = as.numeric(newmeta_x), 
                                                                               y = as.numeric(newmeta_y), stringsAsFactors = F)
                                                                  }
                                                                  else {
                                                                    data.frame(spot = "NA", cell_ratio = "NA", celltype = "NA", 
                                                                               x = "NA", y = "NA", stringsAsFactors = F)
                                                                  }
                                                                }
  return(newmeta)
}

#Put the seeded dec_celltype function back into the package
environment(.generate_newmeta_doParallel_v2) <- asNamespace("SpaTalk")
assignInNamespace(".generate_newmeta_doParallel", .generate_newmeta_doParallel_v2, ns="SpaTalk")

#Assign seed inside newmeta cell 
.generate_newmeta_cell_v2 <- function (newmeta, st_ndata, sc_ndata, sc_celltype, iter_num, 
                                       n_cores, if_doParallel, rand_seed) {
  newmeta_spotname <- unique(newmeta$spot)
  newmeta_cell <- NULL
  cat(crayon::cyan("Generating single-cell data for each spot", 
                   "\n"))
  if (if_doParallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    newmeta_cell <- foreach::foreach(i = 1:length(newmeta_spotname), 
                                     .combine = "rbind", .packages = "Matrix", .export = ".generate_newmeta_spot") %dopar% 
      {
        set.seed(rand_seed)
        spot_name <- newmeta_spotname[i]
        .generate_newmeta_spot(spot_name, newmeta, st_ndata, 
                               sc_ndata, sc_celltype, iter_num)
      }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
  }
  else {
    for (i in 1:length(newmeta_spotname)) {
      spot_name <- newmeta_spotname[i]
      newmeta_spot <- .generate_newmeta_spot(spot_name, 
                                             newmeta, st_ndata, sc_ndata, sc_celltype, iter_num)
      newmeta_cell <- rbind(newmeta_cell, newmeta_spot)
    }
  }
  newmeta_cell$cell <- paste0("C", 1:nrow(newmeta))
  newmeta_cell <- newmeta_cell[, c(8, 4, 5, 3:1, 7, 6)]
  return(newmeta_cell)
}

#Put the seeded dec_celltype function back into the package
environment(.generate_newmeta_cell_v2) <- asNamespace("SpaTalk")
assignInNamespace(".generate_newmeta_cell", .generate_newmeta_cell_v2, ns="SpaTalk")

#Assign seed inside wrapper function for the generate_newmeta
dec_celltype_v2 <- function (object, sc_data, sc_celltype, min_percent = 0.5, min_nFeatures = 10, 
          if_use_normalize_data = T, if_use_hvg = F, if_retain_other_genes = F, 
          if_doParallel = T, use_n_cores = NULL, iter_num = 1000, method = 1, 
          env = "base", anaconda_path = "~/anaconda3", dec_result = NULL, rand_seed) {
  rand_seed_in <- rand_seed
  if (!is(object, "SpaTalk")) {
    stop("Invalid class for object: must be 'SpaTalk'!")
  }
  if_skip_dec_celltype <- object@para$if_skip_dec_celltype
  if (if_skip_dec_celltype) {
    stop("Do not perform dec_celltype() when providing celltype in createSpaTalk()!")
  }
  if (if_doParallel) {
    if (is.null(use_n_cores)) {
      n_cores <- parallel::detectCores()
      n_cores <- floor(n_cores/2)
    }
    else {
      n_cores <- use_n_cores
    }
    n.threads <- n_cores
    if (n_cores == 1) {
      if_doParallel <- F
      n.threads <- 0
    }
  }
  else {
    n_cores <- 1
    n.threads <- 0
  }
  st_data <- object@data[["rawdata"]]
  st_meta <- object@meta[["rawmeta"]]
  st_dist <- .st_dist(st_meta)
  if (min(st_meta$nFeatures) < min_nFeatures) {
    st_meta[st_meta$nFeatures < min_nFeatures, ]$label <- "less nFeatures"
  }
  st_type <- object@para[["st_type"]]
  if (is(sc_data, "data.frame")) {
    sc_data <- methods::as(as.matrix(sc_data), "dgCMatrix")
  }
  if (is(sc_data, "matrix")) {
    sc_data <- methods::as(sc_data, "dgCMatrix")
  }
  if (!is(sc_data, "dgCMatrix")) {
    stop("st_data must be a data.frame or matrix or dgCMatrix!")
  }
  if (!is.character(sc_celltype)) {
    stop("sc_celltype is not a character!")
  }
  if (ncol(sc_data) != length(sc_celltype)) {
    stop("ncol(sc_data) is not consistent with length(sc_celltype)!")
  }
  if (min_percent >= 1 | min_percent <= 0) {
    stop("Please provide a correct min_percent, ranging from 0 to 1!")
  }
  sc_data <- sc_data[which(rowSums(sc_data) > 0), ]
  if (nrow(sc_data) == 0) {
    stop("No expressed genes in sc_data!")
  }
  colnames(sc_data) <- .rename_chr(colnames(sc_data))
  sc_celltype_new <- .rename_chr(sc_celltype)
  warning_info <- .show_warning(sc_celltype, sc_celltype_new)
  if (!is.null(warning_info)) {
    warning(warning_info)
  }
  sc_celltype <- data.frame(cell = colnames(sc_data), celltype = sc_celltype_new, 
                            stringsAsFactors = F)
  object@data$rawndata <- st_data
  if (if_use_normalize_data == T) {
    st_ndata <- .normalize_data(st_data)
    sc_ndata <- .normalize_data(sc_data)
  }
  else {
    st_ndata <- st_data
    sc_ndata <- sc_data
  }
  sc_ndata_raw <- sc_ndata
  genename <- intersect(rownames(st_data), rownames(sc_data))
  if (length(genename) == 0) {
    stop("No overlapped genes between st_data and sc_data!")
  }
  if (method == 1 & is.null(dec_result)) {
    object@data$rawndata <- st_ndata
    st_ndata <- st_ndata[genename, ]
    sc_ndata <- sc_ndata[genename, ]
    if (st_type == "single-cell") {
      cat(crayon::cyan("Performing Non-negative regression for each cell", 
                       "\n"))
    }
    else {
      cat(crayon::cyan("Performing Non-negative regression for each spot", 
                       "\n"))
    }
    if (if_use_hvg) {
      sc_hvg <- .hvg(sc_ndata, sc_celltype)
      if (length(sc_hvg) > 0) {
        st_ndata <- st_ndata[sc_hvg, ]
        sc_ndata <- sc_ndata[sc_hvg, ]
      }
    }
    if (if_doParallel) {
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
    }
    sc_ref <- .generate_ref(sc_ndata, sc_celltype, if_doParallel)
    if (if_doParallel) {
      doParallel::stopImplicitCluster()
      parallel::stopCluster(cl)
    }
    st_nnlm <- NNLM::nnlm(x = as.matrix(sc_ref), y = as.matrix(st_ndata), 
                          method = "lee", loss = "mkl", n.threads = n.threads)
    st_coef <- t(st_nnlm$coefficients)
  }
  if (method == 2 & is.null(dec_result)) {
    cat(crayon::cyan("Using RCTD to deconvolute", "\n"))
    require(spacexr)
    st_coef <- .run_rctd(st_data, sc_data, st_meta, sc_celltype)
  }
  if (method == 3 & is.null(dec_result)) {
    cat(crayon::cyan("Using Seurat to deconvolute", "\n"))
    require(Seurat)
    st_coef <- .run_seurat(st_data, sc_data, sc_celltype)
  }
  if (method == 4 & is.null(dec_result)) {
    cat(crayon::cyan("Using SPOTlight to deconvolute", "\n"))
    require(SPOTlight)
    st_coef <- .run_spotlight(st_data, sc_data, sc_celltype)
  }
  if (method == 5 & is.null(dec_result)) {
    cat(crayon::cyan("Using deconvSeq to deconvolute", "\n"))
    require(deconvSeq)
    st_coef <- .run_deconvSeq(st_data, sc_data, sc_celltype)
  }
  if (method == 6 & is.null(dec_result)) {
    cat(crayon::cyan("Using stereoscope to deconvolute, please install the stereoscope (python package) first!", 
                     "\n"))
    st_coef <- .run_stereoscope(st_data, sc_data, sc_celltype, 
                                env, anaconda_path)
  }
  if (method == 7 & is.null(dec_result)) {
    cat(crayon::cyan("Using cell2location to deconvolute, please install the cell2location (python package) first!", 
                     "\n"))
    require(anndata)
    require(Seurat)
    require(reticulate)
    require(sceasy)
    st_coef <- .run_cell2location(st_data, st_meta, sc_data, 
                                  sc_celltype, env)
  }
  if (!is.null(dec_result)) {
    if (!is.matrix(dec_result)) {
      stop("Please provide a correct dec_result matrix! See demo_dec_result()!")
    }
    dec_colname <- colnames(dec_result)
    dec_colname <- stringr::str_replace_all(dec_colname, 
                                            pattern = "-", replacement = "_")
    if (!all(dec_colname %in% unique(sc_celltype$celltype))) {
      stop("Celltype name in dec_result must be consistent with the names in scRNA-seq reference!")
    }
    dec_rowname <- rownames(dec_result)
    if (!all(st_meta[, 1] == dec_rowname)) {
      stop("Spot/cell name in dec_result must be consistent with the names in st_meta!")
    }
    st_coef <- dec_result
  }
  coef_name <- colnames(st_coef)
  coef_name <- coef_name[order(coef_name)]
  st_coef <- st_coef[, coef_name]
  object@coef <- st_coef
  st_meta <- cbind(st_meta, .coef_nor(st_coef))
  st_ndata <- st_ndata[genename, ]
  sc_ndata <- sc_ndata[genename, ]
  if (st_type == "single-cell") {
    object@meta$rawmeta <- .determine_celltype(st_meta, min_percent)
  }
  else {
    if (if_doParallel) {
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      newmeta <- SpaTalk:::.generate_newmeta_doParallel(st_meta, 
                                              st_dist, min_percent, rand_seed = rand_seed_in)
      doParallel::stopImplicitCluster()
      parallel::stopCluster(cl)
    }
    else {
      newmeta <- .generate_newmeta(st_meta, st_dist, min_percent)
    }
    newmeta_cell <- SpaTalk:::.generate_newmeta_cell(newmeta, st_ndata, 
                                           sc_ndata, sc_celltype, iter_num, n_cores, if_doParallel, rand_seed = rand_seed_in)
    if (if_retain_other_genes) {
      sc_ndata <- sc_ndata_raw
    }
    else {
      if (if_use_hvg) {
        sc_ndata <- sc_ndata_raw
        sc_ndata <- sc_ndata[genename, ]
      }
    }
    newdata <- sc_ndata[, newmeta_cell$cell_id]
    colnames(newdata) <- newmeta_cell$cell
    object@data$newdata <- methods::as(newdata, Class = "dgCMatrix")
    object@meta$newmeta <- newmeta_cell
    st_meta[st_meta$spot %in% newmeta_cell$spot, ]$celltype <- "sure"
    st_dist <- .st_dist(newmeta_cell)
    object@meta$rawmeta <- st_meta
  }
  cat(crayon::green("***Done***", "\n"))
  object@dist <- st_dist
  object@para$min_percent <- min_percent
  object@para$min_nFeatures <- min_nFeatures
  object@para$if_use_normalize_data <- if_use_normalize_data
  object@para$if_use_hvg <- if_use_hvg
  object@para$iter_num <- iter_num
  return(object)
}

#Put the seeded dec_celltype function back into the package
environment(dec_celltype_v2) <- asNamespace("SpaTalk")
assignInNamespace("dec_celltype", dec_celltype_v2, ns="SpaTalk")

#Seeds need to be placed in the following functions
dec_cci_all_v2 <- function (object, n_neighbor = 10, min_pairs = 5, min_pairs_ratio = 0, 
                            per_num = 1000, pvalue = 0.05, co_exp_ratio = 0.1, if_doParallel = T, 
                            use_n_cores = NULL, rand_seed) {
  if (!is(object, "SpaTalk")) {
    stop("Invalid class for object: must be 'SpaTalk'!")
  }
  if (is.null(use_n_cores)) {
    n_cores <- parallel::detectCores()
    n_cores <- n_cores - 2
  }
  else {
    n_cores <- use_n_cores
  }
  if (if_doParallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
  }
  st_meta <- .get_st_meta(object)
  st_data <- .get_st_data(object)
  celltype_dist <- object@dist
  cellname <- unique(st_meta$celltype)
  celltype_pair <- NULL
  for (i in 1:length(cellname)) {
    d1 <- data.frame(celltype_sender = rep(cellname[i], length(cellname)), 
                     celltype_receiver = cellname, stringsAsFactors = F)
    celltype_pair <- rbind(celltype_pair, d1)
  }
  celltype_pair <- celltype_pair[celltype_pair$celltype_sender != 
                                   celltype_pair$celltype_receiver, ]
  cat(paste0("Note: there are ", length(cellname), " cell types and ", 
             nrow(celltype_pair), " pair-wise cell pairs"), "\n")
  pathways <- object@lr_path$pathways
  ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
  cat(crayon::cyan("Begin to find LR pairs", "\n"))
  if (if_doParallel) {
    all_res <- foreach::foreach(i = 1:nrow(celltype_pair), 
                                .packages = c("Matrix", "reshape2"), .export = c(".get_cellpair", 
                                                                                 ".lr_distance", ".get_tf_res", ".get_score")) %dopar% 
      {
        set.seed(rand_seed)
        celltype_sender <- celltype_pair$celltype_sender[i]
        celltype_receiver <- celltype_pair$celltype_receiver[i]
        cell_pair <- .get_cellpair(celltype_dist, st_meta, 
                                   celltype_sender, celltype_receiver, n_neighbor)
        cell_sender <- st_meta[st_meta$celltype == celltype_sender, 
        ]
        cell_receiver <- st_meta[st_meta$celltype == 
                                   celltype_receiver, ]
        cell_pair_all <- nrow(cell_sender) * nrow(cell_receiver)/2
        if (nrow(cell_pair) > min_pairs) {
          if (nrow(cell_pair) > cell_pair_all * min_pairs_ratio) {
            lrdb <- object@lr_path$lrpairs
            lrdb <- .lr_distance(st_data, cell_pair, 
                                 lrdb, celltype_sender, celltype_receiver, 
                                 per_num, pvalue)
            max_hop <- object@para$max_hop
            receptor_tf <- NULL
            if (nrow(lrdb) > 0) {
              receptor_tf <- .get_tf_res(celltype_sender, 
                                         celltype_receiver, lrdb, ggi_tf, cell_pair, 
                                         st_data, max_hop, co_exp_ratio)
              if (!is.null(receptor_tf)) {
                lrdb <- .get_score(lrdb, receptor_tf)
              }
              else {
                lrdb <- NULL
              }
            }
            if (is.data.frame(lrdb)) {
              if (nrow(lrdb) > 0) {
                list(lrdb = lrdb, receptor_tf = receptor_tf, 
                     cell_pair = cell_pair)
              }
            }
          }
        }
      }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    res_receptor_tf <- NULL
    res_lrpair <- NULL
    res_cellpair <- list()
    m <- 0
    for (i in 1:length(all_res)) {
      all_res1 <- all_res[[i]]
      if (!is.null(all_res1)) {
        m <- m + 1
        res_lrpair <- rbind(res_lrpair, all_res1[[1]])
        res_receptor_tf <- rbind(res_receptor_tf, all_res1[[2]])
        res_cellpair[[m]] <- all_res1[[3]]
        names(res_cellpair)[m] <- paste0(unique(all_res1[[1]]$celltype_sender), 
                                         " -- ", unique(all_res1[[1]]$celltype_receiver))
      }
    }
    if (!is.null(res_lrpair)) {
      object@lrpair <- res_lrpair
    }
    if (!is.null(res_receptor_tf)) {
      object@tf <- res_receptor_tf
    }
    if (length(res_cellpair) > 0) {
      object@cellpair <- res_cellpair
    }
  }
  else {
    for (i in 1:nrow(celltype_pair)) {
      celltype_sender <- celltype_pair$celltype_sender[i]
      celltype_receiver <- celltype_pair$celltype_receiver[i]
      cell_pair <- .get_cellpair(celltype_dist, st_meta, 
                                 celltype_sender, celltype_receiver, n_neighbor)
      cell_sender <- st_meta[st_meta$celltype == celltype_sender, 
      ]
      cell_receiver <- st_meta[st_meta$celltype == celltype_receiver, 
      ]
      cell_pair_all <- nrow(cell_sender) * nrow(cell_receiver)/2
      if (nrow(cell_pair) > min_pairs) {
        if (nrow(cell_pair) > cell_pair_all * min_pairs_ratio) {
          lrdb <- object@lr_path$lrpairs
          lrdb <- .lr_distance(st_data, cell_pair, lrdb, 
                               celltype_sender, celltype_receiver, per_num, 
                               pvalue)
          max_hop <- object@para$max_hop
          receptor_tf <- NULL
          if (nrow(lrdb) > 0) {
            receptor_tf <- .get_tf_res(celltype_sender, 
                                       celltype_receiver, lrdb, ggi_tf, cell_pair, 
                                       st_data, max_hop, co_exp_ratio)
            if (!is.null(receptor_tf)) {
              lrdb <- .get_score(lrdb, receptor_tf)
            }
            else {
              lrdb <- "NA"
            }
          }
          lrpair <- object@lrpair
          if (nrow(lrpair) == 0) {
            if (is.data.frame(lrdb)) {
              object@lrpair <- lrdb
            }
          }
          else {
            if (is.data.frame(lrdb)) {
              lrpair <- rbind(lrpair, lrdb)
              object@lrpair <- unique(lrpair)
            }
          }
          tf <- object@tf
          if (nrow(tf) == 0) {
            if (is.data.frame(receptor_tf)) {
              object@tf <- receptor_tf
            }
          }
          else {
            if (is.data.frame(receptor_tf)) {
              tf <- rbind(tf, receptor_tf)
              object@tf <- unique(tf)
            }
          }
          object@cellpair[[paste0(celltype_sender, " -- ", 
                                  celltype_receiver)]] <- cell_pair
        }
      }
      sym <- crayon::combine_styles("bold", "green")
      cat(sym("***Done***"), paste0(celltype_sender, " -- ", 
                                    celltype_receiver), "\n")
    }
  }
  object@para$min_pairs <- min_pairs
  object@para$min_pairs_ratio <- min_pairs_ratio
  object@para$per_num <- per_num
  object@para$pvalue <- pvalue
  object@para$co_exp_ratio <- co_exp_ratio
  return(object)
} #just needs 1 seed

environment(dec_cci_all_v2) <- asNamespace("SpaTalk")
assignInNamespace("dec_cci_all", dec_cci_all_v2, ns="SpaTalk")


#Seed needs to be placed for the score generation sampling
.get_tf_res_doParallel_v2 <- function (celltype_sender, celltype_receiver, lrdb, ggi_tf, cell_pair, 
                                       st_data, max_hop, co_exp_ratio, rand_seed) {
  receptor_name <- unique(lrdb$receptor)
  receptor_tf <- NULL
  receptor_tf <- foreach::foreach(j = 1:length(receptor_name), 
                                  .packages = "Matrix", .combine = "rbind", .export = c(".generate_ggi_res", 
                                                                                        ".generate_tf_gene_all", ".generate_tf_res", ".random_walk")) %dopar% 
    {
      set.seed(rand_seed)
      ggi_res <- .generate_ggi_res(ggi_tf, cell_pair, receptor_name[j], 
                                   st_data, max_hop, co_exp_ratio)
      if (nrow(ggi_res) > 0) {
        tf_gene_all <- .generate_tf_gene_all(ggi_res, 
                                             max_hop)
        tf_gene_all <- data.frame(gene = names(tf_gene_all), 
                                  hop = tf_gene_all, stringsAsFactors = F)
        tf_gene_all_new <- unique(tf_gene_all)
        tf_gene_all <- tf_gene_all_new$hop
        names(tf_gene_all) <- tf_gene_all_new$gene
        ggi_res$dest_tf_enrich <- "NO"
        if (!is.null(tf_gene_all)) {
          ggi_res[ggi_res$dest %in% names(tf_gene_all), 
          ]$dest_tf_enrich <- "YES"
          receptor_tf_temp <- .generate_tf_res(tf_gene_all, 
                                               celltype_sender, celltype_receiver, receptor_name[j], 
                                               ggi_res)
          receptor_tf_temp$score <- .random_walk(receptor_tf_temp, 
                                                 ggi_res)
          receptor_tf_temp
        }
      }
    }
  return(receptor_tf)
}

environment(.get_tf_res_doParallel_v2) <- asNamespace("SpaTalk")
assignInNamespace(".get_tf_res_doParallel", .get_tf_res_doParallel_v2, ns="SpaTalk")

#Needs to pass seed to seeded .get_tf_res_doParallel
dec_cci_v2 <- function (object, celltype_sender, celltype_receiver, n_neighbor = 10, 
                        min_pairs = 5, min_pairs_ratio = 0, per_num = 1000, pvalue = 0.05, 
                        co_exp_ratio = 0.1, if_doParallel = T, use_n_cores = NULL, rand_seed) {
  rand_seed_in <- rand_seed
  if (!is(object, "SpaTalk")) {
    stop("Invalid class for object: must be 'SpaTalk'!")
  }
  if (is.null(use_n_cores)) {
    n_cores <- parallel::detectCores()
    n_cores <- n_cores - 2
  }
  else {
    n_cores <- use_n_cores
  }
  if (if_doParallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
  }
  st_meta <- .get_st_meta(object)
  st_data <- .get_st_data(object)
  celltype_dist <- object@dist
  if (!celltype_sender %in% st_meta$celltype) {
    stop("Please provide a correct celltype_sender")
  }
  if (!celltype_receiver %in% st_meta$celltype) {
    stop("Please provide a correct celltype_receiver")
  }
  cell_pair <- .get_cellpair(celltype_dist, st_meta, celltype_sender, 
                             celltype_receiver, n_neighbor)
  cell_sender <- st_meta[st_meta$celltype == celltype_sender, 
  ]
  cell_receiver <- st_meta[st_meta$celltype == celltype_receiver, 
  ]
  cell_pair_all <- nrow(cell_sender) * nrow(cell_receiver)/2
  if (nrow(cell_pair) <= min_pairs) {
    stop(paste0("Cell pairs found between ", celltype_sender, 
                " and ", celltype_receiver, " less than min_pairs!"))
  }
  if (nrow(cell_pair) <= cell_pair_all * min_pairs_ratio) {
    stop(paste0("Cell pairs found between ", celltype_sender, 
                " and ", celltype_receiver, " less than min_pairs_ratio!"))
  }
  lrdb <- object@lr_path$lrpairs
  pathways <- object@lr_path$pathways
  ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
  cat(crayon::cyan("Begin to find LR pairs", "\n"))
  if (if_doParallel) {
    lrdb <- .lr_distance_doParallel(st_data, cell_pair, lrdb, 
                                    celltype_sender, celltype_receiver, per_num, pvalue)
  }
  else {
    lrdb <- .lr_distance(st_data, cell_pair, lrdb, celltype_sender, 
                         celltype_receiver, per_num, pvalue)
  }
  max_hop <- object@para$max_hop
  receptor_tf <- NULL
  if (nrow(lrdb) > 0) {
    if (if_doParallel) {
      receptor_tf <- SpaTalk:::.get_tf_res_doParallel(celltype_sender, 
                                            celltype_receiver, lrdb, ggi_tf, cell_pair, st_data, 
                                            max_hop, co_exp_ratio, rand_seed = rand_seed_in)
    }
    else {
      receptor_tf <- .get_tf_res(celltype_sender, celltype_receiver, 
                                 lrdb, ggi_tf, cell_pair, st_data, max_hop, co_exp_ratio)
    }
    if (is.null(receptor_tf)) {
      stop(paste0("No LR pairs found between ", celltype_sender, 
                  " and ", celltype_receiver))
    }
    lrdb <- .get_score(lrdb, receptor_tf)
  }
  else {
    stop(paste0("No LR pairs found between ", celltype_sender, 
                " and ", celltype_receiver))
  }
  lrpair <- object@lrpair
  if (nrow(lrpair) == 0) {
    object@lrpair <- lrdb
  }
  else {
    lrpair <- rbind(lrpair, lrdb)
    object@lrpair <- unique(lrpair)
  }
  tf <- object@tf
  if (nrow(tf) == 0) {
    object@tf <- receptor_tf
  }
  else {
    tf <- rbind(tf, receptor_tf)
    object@tf <- unique(tf)
  }
  object@cellpair[[paste0(celltype_sender, " -- ", celltype_receiver)]] <- cell_pair
  if (if_doParallel) {
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
  }
  object@para$min_pairs <- min_pairs
  object@para$min_pairs_ratio <- min_pairs_ratio
  object@para$per_num <- per_num
  object@para$pvalue <- pvalue
  object@para$co_exp_ratio <- co_exp_ratio
  return(object)
}

environment(dec_cci_v2) <- asNamespace("SpaTalk")
assignInNamespace("dec_cci", dec_cci_v2, ns="SpaTalk")
