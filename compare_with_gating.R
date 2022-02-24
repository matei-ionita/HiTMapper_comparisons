library(tidyverse)
# library(plyr)
source("comparison_utils.R")

design <- "CaCo"
# design <- "longitudinal"
algo_name <- "HiTMapper"

features <- read_csv("features/HiTMapper.csv")
gating <- read_csv("Acute_COVID—AP_Acute COVID—Export Statistics—210714_2018.csv")


#########################################################
############### Run the comparison ######################
#########################################################

gating_feat <- process_gating(gating, design)
algo_feat <- process_algo(features, design)
joint_tall <- join_features(gating_feat, algo_feat, algo_name)
labels <- get_labels(joint_tall, algo_name)


pdf("mapper_vs_gating.pdf", width=12, height=9)
ggplot(joint_tall, aes_string(y="gating", x=algo_name)) +
  facet_wrap(~cell_type, scales="free") +
  geom_point(size=2, aes(color=subject)) +
  geom_text(data=labels, aes(label=label, x=0.6*M, y=1.1*M)) +
  geom_segment(data=labels, aes(x=0,y=0,xend=1.2*M,yend=1.2*M), linetype="dashed") +
  theme_bw()
dev.off()
