library(tidyverse)

feat <- read_csv("features/acuteHiTMapper.csv")


to_compare <- setdiff( names(feat) ,  c("subject", "cohort") )
results <- sapply(to_compare, function(cell) {
  w <- wilcox.test(x = feat %>% filter(cohort == "case") %>% pull(cell),
                   y = feat %>% filter(cohort == "control") %>% pull(cell))
  w$p.value 
})
results <- p.adjust(results, method="fdr")
results[ which(results<0.05) ]

feat_tall <- feat %>%
  select(subject, cohort, all_of(to_compare)) %>%
  pivot_longer(all_of(to_compare), names_to="cell_type", values_to="frequency")

ggplot(feat_tall, aes(x=cohort, y=frequency, fill=cohort)) +
  geom_boxplot() +
  facet_wrap(~cell_type, scales="free_y") +
  theme_bw(base_size = 14)