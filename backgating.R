library(tidyverse)
library(wadeTools)
library(ggnewscale)
source("comparison_utils.R")

algo_name <- "HiTMapper"
load(paste0("features/clustering_", algo_name, ".rda"))

data_dir <- "."
dataset <- "acute"
file_path <- paste0(data_dir, "/data_", dataset, ".rda")
load(file_path)

# Drop unwanted markers
pregating <- c("Live")
unreliable <- c("CD16", "CD194")
keep <- setdiff(colnames(data),union(pregating,unreliable))
data <- data[,keep]

data <- data[which(samples=="ND_BM"),]
clustering <- clustering[which(samples=="ND_BM")]

###############################################
# Get clusters based on labels
###############################################
cd45 <- seq_along(clustering)
tcells <- grep("Tcells", clustering)
cd4 <- grep("CD4", clustering)
cd4_naive <- grep("CD4Tcells_Naive", clustering)
treg <- grep("Treg", clustering)
tfh <- grep("Tfh", clustering)

cd8 <- grep("CD8", clustering)
cd8_naive <- grep("CD8Tcells_Naive", clustering)

gd <- grep("gdT", clustering)

bcells <- grep("Bcells", clustering)
pb <- grep("PB", clustering)

baso <- grep("Baso", clustering)
eosi <- grep("Eosi", clustering)
neut <- grep("Neut", clustering)

mono <- grep("Mono", clustering)
mono_nc <- grep("NCMono", clustering)
mono_cl <- setdiff(mono, mono_nc)

pdc <- grep("pDC", clustering)



###############################################
# plot
###############################################
plot_path <- "figures/backgating/ND_BM/"

save_backgate(data, params=c("CD3", "CD19"), background=cd45, cluster=tcells,
              path_base=paste0(plot_path, "backgate_tcell_", algo_name))
save_backgate(data, params=c("CD4", "CD8a"), background=tcells, cluster=cd4,
              path_base=paste0(plot_path, "backgate_cd4_", algo_name))
save_backgate(data, params=c("CD4", "CD8a"), background=tcells, cluster=cd8,
              path_base=paste0(plot_path, "backgate_cd8_", algo_name))

save_backgate(data, params=c("CD4", "TCRgd"), background=tcells, cluster=gd,
              path_base=paste0(plot_path, "backgate_gd_", algo_name))

save_backgate(data, params=c("CD185", "CD4"), background=cd4, cluster=tfh,
              path_base=paste0(plot_path, "backgate_tfh_", algo_name))
save_backgate(data, params=c("CD25", "CD127"), background=cd4, cluster=treg,
              path_base=paste0(plot_path, "backgate_treg_", algo_name))


save_backgate(data, params=c("CD3", "CD19"), background=cd45, cluster=bcells,
              path_base=paste0(plot_path, "backgate_bcell_", algo_name))
save_backgate(data, params=c("CD27", "CD38"), background=bcells, cluster=pb,
              path_base=paste0(plot_path, "backgate_pb_", algo_name))

save_backgate(data, params=c("CD66b", "CD294"), background=cd45, cluster=neut,
              path_base=paste0(plot_path, "backgate_neut_", algo_name))
save_backgate(data, params=c("CD66b", "CD294"), background=cd45, cluster=eosi,
              path_base=paste0(plot_path, "backgate_eosi_", algo_name))
save_backgate(data, params=c("CD66b", "CD294"), background=cd45, cluster=baso,
              path_base=paste0(plot_path, "backgate_baso_", algo_name))
save_backgate(data, params=c("CD123", "CD294"), background=cd45, cluster=baso,
              path_base=paste0(plot_path, "backgate_baso_", algo_name))


save_backgate(data, params=c("CD123", "HLA-DR"), background=cd45, cluster=pdc,
              path_base=paste0(plot_path, "backgate_pdc_", algo_name))

save_backgate(data, params=c("CD11c", "CD66b"), background=cd45, cluster=mono,
              path_base=paste0(plot_path, "backgate_mono_", algo_name))
save_backgate(data, params=c("CD11c", "CD14"), background=mono, cluster=mono_cl,
              path_base=paste0(plot_path, "backgate_mono_cl_", algo_name))
save_backgate(data, params=c("CD11c", "CD14"), background=mono, cluster=mono_nc,
              path_base=paste0(plot_path, "backgate_mono_nc_", algo_name))




backgate(data, c("CD4", "CD8a"), background = tcells, cluster=cd4)

ff <- flowFrame(exprs = data[,c("CD45","CD66b")])
ggflow(ff, params=c("CD45", "CD66b"), trans_fl = "asinh")

