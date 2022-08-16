library(HiTMapper)
library(tidyverse)

`%+%` <- paste0
cohort <- "acute"
data_dir <- "data/"
load(data_dir %+% "data_" %+% cohort %+% ".rda")

# Drop unwanted markers
pregating <- c("Live")
unreliable <- c("CD16", "CD194")
keep <- setdiff(colnames(data), c(unreliable, pregating))

####################################################
################ Run algorithm #####################
####################################################
start_time <- Sys.time()
set.seed(127)
mapper <- HiTMapper(data[,keep], total_nodes=1000, 
                    grid_size=c(9,8),
                    outlier_cutoff = 100,
                    resolution=8)
end_time <- Sys.time()
print(end_time-start_time)

defs <- read_csv("maxpar_phenos_functional.csv")
mapper <- extract_features(data[,keep], samples, mapper)
mapper <- label_communities(mapper, defs)

save(mapper, file="mapper_" %+% cohort %+% ".rda")

clustering <- mapper$community_mapping
save(clustering, file = "features/clustering_" %+% cohort %+% "HiTMapper.rda")

####################################################
####### Standardize features and write to file #####
####################################################

cell_types <- levels(mapper$community) %>%
  as.character() %>%
  str_split(fixed(".")) %>%
  sapply(function(x) x[1]) %>%
  unique()

feat <- lapply(cell_types , function(un) {
  sel_cols <- grep(un, colnames(mapper$features))
  apply(mapper$features[,sel_cols,drop=FALSE], 1, sum)
}) %>% do.call(what=cbind)
colnames(feat) <- cell_types
feat <- as_tibble(feat)

feat$subject <- row.names(mapper$features)
feat$cohort <- row.names(mapper$features) %>%
  sapply(function(name) {
    if(substr(name,1,2) %in% c("HD","ND"))
      return("control")
    else
      return("case")
  })
feat <- feat %>% relocate(subject, cohort)
write_csv(feat, "features/" %+% cohort %+% "HiTMapper.csv")

