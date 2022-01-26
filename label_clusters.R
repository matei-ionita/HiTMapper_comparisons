
###############################################
# assuming that cluster_mapping is a factor
# and sample_mapping is a character vector...
# fix later?
###############################################
extract_features_from_clusters <- function(cluster_mapping, sample_mapping, data) {
  # features
  base <- unname(table(sample_mapping)) %>% as.numeric()
  tab <- table(cluster_mapping, sample_mapping) %>% as.matrix()
  feat <- apply(tab, 1, function(row) row/base)
  if (length(base) == 1)
    feat <- as.matrix(feat) %>% t()

  # cluster centroids
  comms <- levels(cluster_mapping)
  comm_events <- lapply(comms, function(comm) which(cluster_mapping==comm))
  centroids <- sapply(comm_events, function(events) {
    apply(data[events,,drop=FALSE],2,median)
  }) %>% t()

  return(list(features=feat, centroids=centroids))
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




####################################################
# Extract metadata from sample names
####################################################
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



