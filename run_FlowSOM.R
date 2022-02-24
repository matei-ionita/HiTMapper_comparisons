library(FlowSOM)
library(tidyverse)
source("label_clusters.R")


data_dir <- "."
dataset <- "acute"
file_path <- paste0(data_dir, "/data_", dataset, ".rda")
load(file_path)

# Drop unwanted markers
pregating <- c("Live")
unreliable <- c("CD16", "CD194")
keep <- setdiff(colnames(data),union(pregating,unreliable))
data <- data[,keep]


# Just for a quick trial run, subsample the data
set.seed(4001)
sel <- sample(nrow(data), 5e5)
dat <- data[sel,]
sam <- samples[sel]


####################################################
############### Run the algorithm ##################
####################################################
start_time <- Sys.time()
set.seed(4003)
# For the final run (with all data), try something
# like 25x25 grid, with nClus=50 or 60?
fsom <- FlowSOM(dat, colsToUse = colnames(dat),
                xdim=25, ydim=25, nClus=50)
cluster_mapping <- GetMetaclusters(fsom)
end_time <- Sys.time()
print(end_time-start_time)

####################################################
# Post-processing: feature extraction and labeling #
####################################################
defs <- read_csv("maxpar_phenos_functional.csv")

post <- extract_features_from_clusters(cluster_mapping, sam, dat)
features <- data.frame(post$features)
labels <- label_clusters(post$centroids, defs)
names(features) <- labels

meta <- get_meta(row.names(features), design="CaCo")
features <- cbind(meta,features)

write_csv(features, "features/FlowSOM.csv")



