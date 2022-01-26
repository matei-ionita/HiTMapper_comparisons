
process_gating <- function(gating, design) {
  keep <- which(grepl("percent_total", names(gating)))
  gating_feat <- gating[,keep] / 100
  gating_feat <- apply(gating_feat, 2, function(x) x/gating_feat$`CD45+|percent_total`) %>%
    as_tibble()
  names(gating_feat) <- str_remove_all(names(gating_feat), fixed("|percent_total"))
  gating_feat$`CD8 mem` <- gating_feat[,which(grepl("CD8 TEM", names(gating_feat)) | grepl("CD8 TCM", names(gating_feat)))] %>%
    apply(1, sum)
  gating_feat$`CD4 mem` <- gating_feat[,which(grepl("CD4 TEM", names(gating_feat)) | grepl("CD4 TCM", names(gating_feat)))] %>%
    apply(1, sum)
  gating_feat$`B cells mem` <- gating_feat$`Total Memory B cells`
  gating_feat$`B cells Naive` <- gating_feat$`Naive IgD+`
  gating_feat$`B cells CD27neg` <- gating_feat$`CD27-B`

  meta_gating <- get_sample_time(gating$file %>% str_remove("_Normalized.fcs"),
                                 design=design)
  gating_feat <- cbind(meta_gating, gating_feat)

  return(gating_feat)
}

process_algo <- function(features, design) {
  cols <- names(features)

  bcell <- which(grepl("Bcells", cols))
  bcell_naive <- which(grepl("Bcells_Naive", cols))
  bcell_27neg <- which(grepl("Bcells_Naive", cols) |
                         grepl("CD27negIgDneg", cols))
  bcell_mem <- which(grepl("switch", cols))
  bcell_pb <- which(grepl("Bcells_PB", cols))

  tcell <- which(grepl("Tcells", cols))

  tcell_cd8 <- which(grepl("CD8", cols))
  tcell_cd8_naive <- which(grepl("CD8Tcells_Naive", cols))
  tcell_cd8_mem <- setdiff(tcell_cd8, tcell_cd8_naive)
  tcell_cd8_mem_act <- which(grepl("CD8Tcells_Mem", cols) & grepl("act", cols))

  tcell_cd4 <- which(grepl("CD4", cols))
  tcell_cd4_naive <- which(grepl("CD4Tcells_Naive", cols))
  tcell_cd4_mem <- setdiff(tcell_cd4, tcell_cd4_naive)
  tcell_cd4_mem_act <- which(grepl("CD4Tcells", cols) & grepl("act", cols))

  tfh <- which(grepl("Tfh", cols))
  th1 <- which(grepl("Th1", cols) & !grepl("Th17", cols))
  th2 <- which(grepl("Th2", cols))
  # th17 <- which(grepl("Th17", cols) & !grepl("Th1", cols))
  th17 <- which(grepl("Th17", cols))

  tcell_gd <- which(grepl("gd", cols))
  nk <- which(grepl("Nkcells", cols))
  lymph <- union_all(bcell, tcell, nk)

  pdc <- which(grepl("pDC", cols))
  mdc <- which(grepl("mDC", cols))
  mono <- which(grepl("Monocyte", cols))

  eos <- which(grepl("Eosinophils", cols))
  bas <- which(grepl("Basophils", cols))
  neu <- which(grepl("Neutrophil", cols))

  algo_feat <- tibble(Lymphocytes = apply(features[,lymph,drop=FALSE], 1, sum),
                        `B cells` = apply(features[,bcell,drop=FALSE], 1, sum),
                        `B cells mem` = apply(features[,bcell_mem,drop=FALSE], 1, sum),
                        `B cells Naive` = apply(features[,bcell_naive,drop=FALSE], 1, sum),
                        `B cells CD27neg` = apply(features[,bcell_27neg,drop=FALSE], 1, sum),
                        Plasmablasts = apply(features[,bcell_pb,drop=FALSE], 1, sum),
                        `T cells` = apply(features[,tcell,drop=FALSE], 1, sum),
                        `CD8+ T` = apply(features[,tcell_cd8,drop=FALSE], 1, sum),
                        `CD8 Naive` = apply(features[,tcell_cd8_naive,drop=FALSE], 1, sum),
                        `CD8 mem` = apply(features[,tcell_cd8_mem,drop=FALSE], 1, sum),
                        `activated nnCD8 T` = apply(features[,tcell_cd8_mem_act,drop=FALSE], 1, sum),
                        `CD4+ T` = apply(features[,tcell_cd4,drop=FALSE], 1, sum),
                        `CD4 Naive` = apply(features[,tcell_cd4_naive,drop=FALSE], 1, sum),
                        `CD4 mem` = apply(features[,tcell_cd4_mem,drop=FALSE], 1, sum),
                        `activated nnCD4 T` = apply(features[,tcell_cd4_mem_act,drop=FALSE], 1, sum),
                        `NNCD4 T CXCR5+` = apply(features[,tfh,drop=FALSE], 1, sum),
                        `Th1` = apply(features[,th1,drop=FALSE], 1, sum),
                        `Th2` = apply(features[,th2,drop=FALSE], 1, sum),
                        `Th17` = apply(features[,th17,drop=FALSE], 1, sum),
                        `gd T cells` = apply(features[,tcell_gd,drop=FALSE], 1, sum),
                        `Total Monocytes` = apply(features[,mono,drop=FALSE], 1, sum),
                        `pDC` = apply(features[,pdc,drop=FALSE], 1, sum),
                        `mDC` = apply(features[,mdc,drop=FALSE], 1, sum),
                        `Total NK cells` = apply(features[,nk,drop=FALSE], 1, sum),
                        Basophils = apply(features[,bas,drop=FALSE], 1, sum),
                        Eosinophils = apply(features[,eos,drop=FALSE], 1, sum),
                        Neutrophils = apply(features[,neu,drop=FALSE], 1, sum)
  )


  if (design=="CaCo")
    meta <- data.frame(subject = features$subject,
                       cohort = features$cohort)
  else
    meta <- data.frame(subject = features$subject,
                       time = features$time)
  algo_feat <- cbind(meta, algo_feat)
  return(algo_feat)
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


get_sample_time <- function(sample_names, design) {
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


