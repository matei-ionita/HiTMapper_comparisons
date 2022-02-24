library(HiTMapper)
library(tidyverse)

data_dir <- "."
dataset <- "acute"
file_path <- paste0(data_dir, "/data_", dataset, ".rda")
load(file_path)

# Drop unwanted markers
pregating <- c("Live")
unreliable <- c("CD16", "CD194")
keep <- setdiff(colnames(data),union(pregating,unreliable))
data <- data[,keep]
# data[,c("CD38","HLA-DR","CD185")] <- 2*data[,c("CD38","HLA-DR","CD185")]

####################################################
################ Run algorithm #####################
####################################################
start_time <- Sys.time()
set.seed(127)
mapper <- HiTMapper(data, total_nodes=1000, outlier_cutoff = 100)
end_time <- Sys.time()
print(end_time-start_time)

defs <- read_csv("maxpar_phenos_functional.csv")
mapper <- extract_features(data, samples, mapper)
mapper <- label_communities(mapper, defs)



save(mapper, file="mapper.rda")

clustering <- mapper$community_mapping
save(clustering, file = "features/clustering_HiTMapper.rda")

####################################################
####### Standardize features and write to file #####
####################################################
features <- data.frame(mapper$features)
features$subject <- row.names(mapper$features)
features$cohort <- row.names(mapper$features) %>%
  sapply(function(name) {
    if(substr(name,1,2) %in% c("HD","ND"))
      return("control")
    else
      return("case")
  })
row.names(features) <- NULL
write_csv(features, "features/HiTMapper.csv")

