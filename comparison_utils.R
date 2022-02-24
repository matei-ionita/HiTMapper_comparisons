
process_gating <- function(gating, design) {
  keep <- which(grepl("percent_total", names(gating)))
  gating_feat <- gating[,keep] / 100
  gating_feat <- apply(gating_feat, 2, function(x) x/gating_feat$`CD45+|percent_total`) %>%
    as_tibble()
  names(gating_feat) <- str_remove_all(names(gating_feat), fixed("|percent_total"))

  cd8mem <- which(grepl("CD8 TEM", names(gating_feat)) | grepl("CD8 TCM", names(gating_feat)))
  gating_feat$`CD8 mem` <- gating_feat[,cd8mem] %>% apply(1, sum)

  cd4mem <- which(grepl("CD4 TEM", names(gating_feat)) | grepl("CD4 TCM", names(gating_feat)))
  gating_feat$`CD4 mem` <- gating_feat[,cd4mem] %>% apply(1, sum)

  gating_feat$`B cells mem` <- gating_feat$`Total Memory B cells`
  gating_feat$`B cells Naive` <- gating_feat$`Naive IgD+`
  gating_feat$`B cells CD27neg` <- gating_feat$`CD27-B`

  file_names <- gating$file %>% str_remove("_Normalized.fcs")
  meta_gating <- get_meta(file_names,design=design)
  gating_feat <- cbind(meta_gating, gating_feat)

  return(gating_feat)
}

process_algo <- function(features, design) {
  cols <- names(features)
  col_groups <- aggregate_cols(cols)

  algo_feat <- sapply(col_groups, function(group) {
    apply(features[,group,drop=FALSE], 1, sum)
  })

  if (design=="CaCo")
    meta <- data.frame(subject = features$subject,
                       cohort = features$cohort)
  else
    meta <- data.frame(subject = features$subject,
                       time = features$time)

  algo_feat <- cbind(meta, algo_feat)
  return(algo_feat)
}


aggregate_cols <- function(cols) {
  cols <- cols %>% str_split(pattern=fixed(".")) %>%
    sapply(function(x) x[1])

  groups <- list()
  groups[["B cells"]] <- which(grepl("Bcells", cols))
  groups[["B cells Naive"]] <- which(grepl("Bcells_Naive", cols))
  groups[["B cells CD27neg"]] <- which(grepl("Bcells_Naive", cols) |
                         grepl("CD27negIgDneg", cols))
  groups[["B cells mem"]] <- which(grepl("switch", cols))
  groups[["Plasmablasts"]] <- which(grepl("Bcells_PB", cols))

  groups[["T cells"]] <- which(grepl("Tcells", cols))

  groups[["CD8+ T"]] <- which(grepl("CD8", cols))
  groups[["CD8 Naive"]] <- which(grepl("CD8Tcells_Naive", cols))
  groups[["CD8 mem"]] <- which(grepl("CD8", cols) & !grepl("Naive", cols))
  groups[["activated nnCD8 T"]] <- which(grepl("CD8Tcells_Mem", cols) & grepl("act", cols))

  groups[["CD4+ T"]] <- which(grepl("CD4", cols))
  groups[["CD4 Naive"]] <- which(grepl("CD4Tcells_Naive", cols))
  groups[["CD4 mem"]] <- which(grepl("CD4", cols) & !grepl("Naive", cols))
  groups[["activated nnCD4 T"]] <- which(grepl("CD4Tcells", cols) & grepl("act", cols))

  groups[["NNCD4 T CXCR5+"]] <- which(grepl("Tfh", cols))
  groups[["Th1"]] <- which(grepl("Th1", cols) & !grepl("Th17", cols))
  groups[["Th2"]] <- which(grepl("Th2", cols))
  # th17 <- which(grepl("Th17", cols) & !grepl("Th1", cols))
  groups[["Th17"]] <- which(grepl("Th17", cols))

  groups[["gd T cells"]] <- which(grepl("gd", cols))
  groups[["Total NK cells"]] <- which(grepl("Nkcells", cols))
  groups[["Lymphocytes"]] <- union(union(groups[["B cells"]], groups[["T cells"]]),
                                   groups[["Total NK cells"]])

  groups[["pDC"]] <- which(grepl("pDC", cols))
  groups[["mDC"]] <- which(grepl("mDC", cols))
  groups[["Total Monocytes"]] <- which(grepl("Monocyte", cols))

  groups[["Eosinophils"]] <- which(grepl("Eosinophils", cols))
  groups[["Basophils"]] <- which(grepl("Basophils", cols))
  groups[["Neutrophils"]] <- which(grepl("Neutrophil", cols))
  return(groups)
}


join_features <- function(gating_feat, algo_feat, algo_name) {
  keep <- intersect(names(gating_feat), names(algo_feat))

  gating_tall <- select_and_pivot(gating_feat, keep, "gating")
  algo_tall <- select_and_pivot(algo_feat, keep, algo_name)

  joint <- rbind(gating_tall, algo_tall) %>%
    pivot_wider(names_from = "method", values_from = "fraction")

  return(joint)
}

select_and_pivot <- function(feat, keep, method) {
  phenos <- setdiff(keep, c("subject", "time", "cohort"))
  tall <- feat %>%
    select(all_of(keep)) %>%
    pivot_longer(all_of(phenos), names_to = "cell_type",
                 values_to = "fraction")
  tall$method <- method
  return(tall)
}


get_meta <- function(sample_names, design) {
  if (design == "CaCo") {
    meta <- tibble(subject = sample_names,
                   cohort = sapply(sample_names, get_cohort))
    return(meta)
  }

  meta <- tibble(subject = sample_names %>%
                   str_split("[._]+") %>%
                   sapply(get_sample),
                 time = sample_names %>%
                   str_split("[._]+") %>%
                   sapply(get_time))

  return(meta)
}

get_cohort <- function(name) {
  if(substr(name,1,2) %in% c("HD","ND"))
    return("control")
  else
    return("case")
}

get_sample <- function(x) {
  if (substr(x[[1]],1,1) == "H")
    return("HD")
  return(x[[1]])
}

get_time <- function(x) {
  if (length(x) == 2)
    return(x[[2]])
  if (substr(x[[1]],1,1) == "H")
    return(x[[1]] %>% str_remove("HD") %>% str_remove("H"))
  return(paste0(x[[2]], "_", x[[3]]))
}


get_labels <- function(joint_tall, algo_name) {
  phenos <- unique(joint_tall$cell_type)
  rows <- lapply(phenos, function(pheno) which(joint_tall$cell_type==pheno))
  names(rows) <- phenos

  gating <- joint_tall$gating
  algo <- pull(joint_tall, algo_name)

  labels <- tibble(cell_type=phenos,
                   M=sapply(phenos, function(pheno) {
                     max(c(gating[rows[[pheno]]],
                           algo[rows[[pheno]]]))
                   }),
                   label=sapply(phenos, function(pheno) {
                     x <- algo[rows[[pheno]]]
                     y <- gating[rows[[pheno]]]
                     paste("rho =", round(cor(x, y, method="spearman"),3))
                   }))
}




########################################################
# BACKGATING
########################################################


backgate <- function(data, params, background, cluster) {
  background <- setdiff(background, cluster)

  js <- js_divergence(x=data[background,params],
                      y=data[cluster,params]) %>%
    round(digits=3)

  x_max <- max(7, data[c(background,cluster),params[1]])
  y_max <- max(7, data[c(background,cluster),params[2]])

  ff <- flowFrame(exprs=data[background,params])
  p <- ggflow(ff, params=params, trans_fl = "asinh",
              xlim=c(0,x_max),
              ylim=c(0,y_max))

  df <- data.frame(data[cluster,params], check.names=FALSE)

  p+new_scale_fill()+
    geom_bin2d(data=df, mapping=aes(.data[[params[1]]],.data[[params[2]]]),
               bins=750,na.rm=TRUE) +
    viridis::scale_fill_viridis(option = "rocket") +
    annotate("text", label=paste("JS divergence:", js),
             x=x_max/2, y=y_max-0.5, size=10)
}

save_backgate <- function(data, params, background, cluster,
                          path_base) {
  p <-  backgate(data, params, background, cluster)

  png(paste0(path_base, "_", params[1], "_", params[2], ".png"),
      width=800, height=700)
  plot(p)
  dev.off()
}


standardize <- function(x, tol=1e-12) {
  m <- mean(x)
  s <- sd(x)

  if(s < tol)
    return(x-m)
  return((x-m)/s)
}


js_divergence <- function(x,y,bins=256) {
  lims <- list(c(min(c(x[,1],y[,1])),max(c(x[,1],y[,1]))),
               c(min(c(x[,2],y[,2])),max(c(x[,2],y[,2]))))
  bw <- vapply(lims, function(lim) (lim[2]-lim[1])/bins, numeric(1))/2
  p <- KernSmooth::bkde2D(x,bandwidth=bw,
                          gridsize=c(bins,bins),
                          range.x = lims)$fhat
  p <- p/sum(p)
  q <- KernSmooth::bkde2D(y,bandwidth=bw,
                          gridsize=c(bins,bins),
                          range.x = lims)$fhat
  q <- q/sum(q)
  m <- 0.5*(p+q)
  KL1 <- p*log(p/m,base=2)
  KL2 <- q*log(q/m,base=2)

  js <- 0.5*sum(KL1[which(p>0)], KL2[which(q>0)])
}




