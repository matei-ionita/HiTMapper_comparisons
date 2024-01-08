get_centroids <- function(clustering, data) {

  comms <- levels(clustering)
  comm_events <- lapply(comms, function(comm) which(clustering==comm))
  centroids <- sapply(comm_events, function(events) {
    apply(data[events,,drop=FALSE],2,median)
  }) %>% t()
  
  return(centroids)
}


label_clusters <- function(centroids, defs) {
  modality <- apply(centroids, 2, get_modality)
  matches <- match_defs(defs, modality)
  phenos <- c(defs$Phenotype, "Other")
  labels <- phenos[matches]
  labels <- make.unique(labels)
  
  return(labels)
}

##########################################
# Divisive hierarchical clustering (diana)
# of the cluster centroids into "hi" and
# "lo" modalities.
# Independently for all markers.
##########################################
get_modality <- function(marker)
{
  diana <- cluster::diana(marker)
  diana_cut <- cutree(as.hclust(diana), k = 2)
  
  if (mean(marker[which(diana_cut == 1)]) < mean(marker[which(diana_cut == 2)])) {
    mod <- c("lo", "hi")
  }
  else {
    mod <- c("hi", "lo")
  }
  return(mod[diana_cut])
}


#########################################
# For each phenotype, check for clusters
# whose modality matches the definition,
# and label them.
# IMPORTANT: this process is sequential,
# and matches down the list overwrite
# previous ones. Design the phenotype
# spreadsheet with this in mind.
#########################################
match_defs <- function(defs, modality) {
  ind <- integer(nrow(modality))
  
  for (i in seq(nrow(defs))) {
    markers <- names(defs)[which(defs[i,] %in% c("hi", "lo"))]
    match <- apply(modality[,markers,drop=FALSE], 1, function(x) {
      def <- as.character(defs[i,markers])
      chx <- as.character(x)
      all.equal(chx,def)
    })
    sel <- which(match == "TRUE")
    ind[sel] <- i
  }
  
  ind[which(ind==0)] <- nrow(defs)+1
  return(ind)
}


get_intermediate_populations <- function(feat) {
  feat %>% 
    mutate(`T cells CD4` = `T cells CD4 Naive` + `T cells CD4 Tfh` +
             `T cells CD4 Treg` + `T cells CD4 Memory Th1` +
             `T cells CD4 Memory not Th1`,
           `T cells CD8` = `T cells CD8 Naive` + `T cells CD8 Memory`) %>%
    mutate(`T cells` = `T cells CD4` + `T cells CD8` + `T cells MAIT NKT` +
             `T cells gd`,
           `Monocytes` = `Monocytes Classical` + `Monocytes Nonclassical`)
}


rename_features <- function(gating) {
  gating %>%
    mutate(`T cells CD4` = `CD4+ T`,
           `T cells CD8` = `CD8+ T`,
           `T cells CD4 Naive` = `CD4 Naive`,
           `T cells CD8 Naive` = `CD8 Naive`,
           `T cells CD4 Tfh` = `NNCD4 T CXCR5+`,
           `T cells CD4 Memory Th1` = `Th1`,
           `T cells gd` = `gd T cells`)
}


