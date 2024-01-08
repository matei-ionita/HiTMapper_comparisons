library(tidyverse)
library(HiTMapper)
library(FlowSOM)
library(FastPG)
library(uwot)

## modify paths to fit your local dir structure
# dir_data <- "data/vaccine_covid_cytof"
dir_data <- "data/acute_covid_flow"
dir_script <- "scripts"
dir_result <- "results"
dir_figure <- "plots/acute_flow"
dir_phenos <- "phenos"

source(paste0(dir_script, "/utils.R"))


load(paste0(dir_data, "/data.rda"))
phenos <- read_csv(paste0(dir_phenos, "/flow_phenos_lineage.csv"))
metadata <- read_csv(paste0(dir_data, "/OMIQ_metadata-Acute_COVID_flow.csv")) %>%
  select(-OmiqID)
  # mutate(Vaccine = str_remove(Vaccine, ","),
  #        Control = grepl("H", Filename))

#### HiTMapper #####
set.seed(0)
system.time(mapper <- HiTMapper(data, defs=phenos, resolution = 4))
save(mapper, file=paste0(dir_result, "/mapper_acute_flow.rda"))


#### FlowSOM #####
## optional, compare to FlowSOM results
set.seed(0)
system.time(fsom <- FlowSOM(data, nClus=12))
save(fsom, file=paste0(dir_result, "/fsom_acute_flow.rda"))
flowSOM_clustering <- GetMetaclusters(fsom)

fsom_centroids <- get_centroids(flowSOM_clustering, data)
fsom_labels <- label_clusters(fsom_centroids, phenos)
fsom_clustering <- as.factor(fsom_labels[flowSOM_clustering])


#### FastPG #####
## optional, compare to FastPG results
## warning: much slower than HiTMapper or FlowSOM
set.seed(0)
system.time(fastpg <- fastCluster(data, verbose = TRUE))

pg_centroids <- get_centroids(as.factor(fastpg$communities), data)
pg_labels <- label_clusters(pg_centroids, phenos)
pg_clustering <- as.factor(pg_labels[fastpg$communities + 1])

#### umap for visualization #####

set.seed(0)
sel_umap <- sample(nrow(data), 2e5)
um <- umap(data[sel_umap,], verbose=TRUE)
# save(um, file=paste0(dir_result, "/umap_acute_flow.rda"))
# load(file=paste0(dir_result, "/umap_acute_flow.rda"))

df <- as_tibble(data[sel_umap,]) %>%
  mutate(UMAP1 = um[,1], 
         UMAP2 = um[,2],
         cluster_h = mapper$clustering[sel_umap],
         cluster_f = fsom_clustering[sel_umap],
         cluster_p = pg_clustering[sel_umap],
         Filename = samples[sel_umap]) %>%
  inner_join(metadata)

ggplot(df, aes(x=UMAP1, y=UMAP2, color=cluster_h)) + 
  geom_point(shape=1, size=0.7, alpha=0.5) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, shape=19))) +  
  theme_bw()
ggsave(filename = paste0(dir_figure, "/umap/clusters.png"), width=8, height=5)

for (m in colnames(data)) {
  g <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(aes(color=.data[[m]]), shape=1, size=0.5, alpha=0.5) +
    scale_color_gradient2(low="blue", mid="black", high="red") +
    theme_bw(base_size=12)
  ggsave(paste0(dir_figure, "/umap/channel_", m, ".png"),
         plot = g, width = 6, height=5)
}



