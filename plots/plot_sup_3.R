library(tidyverse)
library(patchwork)

samples <- c(
  "dnm1l_barcode01", "gridION_Q20", "gridION",
  "longbench_bulk_H146",
  "longbench_sc_ont", "longbench_sc_pb", "lr-split-seq",
  "promethION", "rage-seq", "scmixology2", "visium"
)

df <- samples %>%
  set_names() %>%
  map(~ read_csv(file.path("data/fdd_stats", paste0(.x, ".csv")), show_col_types = FALSE)) %>%
  list_rbind(names_to = "dataset") %>%
  mutate(
    ratio = new_nodes / valid_nodes,
    dataset = factor(dataset, levels = samples)
  )

df.subset <- df %>%
  filter(valid_nodes > 25)
# slice_sample(n = 100000)

p_all <- ggplot(df.subset, aes(x = ratio, colour = dataset)) +
  geom_density() +
  coord_cartesian(ylim = c(0, 25)) +
  geom_vline(xintercept = 0.25) +
  labs(title = "All datasets", x = "structural divergence ratio", y = "density") +
  guides(colour = guide_legend(nrow = 2)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

p_each <- ggplot(df.subset, aes(x = ratio)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5, bins = 100) +
  geom_vline(xintercept = 0.25) +
  facet_wrap(~dataset, ncol = 4, scales = "free_y") +
  coord_cartesian(xlim = c(0, 0.6)) +
  labs(x = "structural divergence ratio", y = "density") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

plot <- p_all / p_each + plot_layout(heights = c(1.5, 3.5))

ggsave("results/figures/sup_fig_3.pdf", plot, width = 210, height = 220, units = "mm")
cat("Done\n")
