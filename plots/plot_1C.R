# Builds Figure 1C
library(tidyverse)

df.clustering_test_dataset <- read_tsv("visium_header_stats.tsv", show_col_types = FALSE) %>%
  mutate(
    ratio = new_nodes / valid_nodes,
    predicted_type = factor(ifelse((ratio < 0.25) & (valid_nodes > 25), "DUP", "RND")),
    real_type = factor(type, levels = c("DUP", "RND"))
  ) %>%
  select(-type)

source("plots/subplots/fdd_ratio.R")

# compute precision, recall
confusion <- df.clustering_test_dataset %>%
  summarise(
    tp = sum(real_type == "DUP" & predicted_type == "DUP"),
    fp = sum(real_type == "RND" & predicted_type == "DUP"),
    fn = sum(real_type == "DUP" & predicted_type == "RND"),
    tn = sum(real_type == "RND" & predicted_type == "RND")
  )
precision <- with(confusion, tp / (tp + fp))
recall <- with(confusion, tp / (tp + fn))
f1 <- 2 * (precision * recall) / (precision + recall)
cat("Precision: ", precision)
cat("Recall: ", recall)
cat("F1: ", f1)

# Dimensions match the manuscript column width.
ggsave("results/figures/fig_1c.png", fig.fdd_ratio, width = 90, height = 40, units = "mm", dpi = 300)
ggsave("results/figures/fig_1c.svg", fig.fdd_ratio, width = 90, height = 40, units = "mm")
