# ==== DATA ====

df.clustering_test_dataset %>%
  count(real_type, predicted_type) %>%
  write_csv("results/tables/fdd_ratio.csv")


# ==== PLOT ====

fig.fdd_ratio <- df.clustering_test_dataset %>%
  filter(valid_nodes > 25 & ratio < 0.5) %>%
  ggplot(aes(x = ratio, fill = real_type)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.5,
    bins = 100
  ) +
  # draw threshold cutoff
  geom_vline(xintercept = 0.25) +
  annotate("text", label = "threshold", x = 0.3, y = 17) +
  # scale_y_continuous(labels = NULL) +
  scale_x_continuous(limits = c(0, 0.5)) +
  scale_fill_manual(
    values = c("DUP" = dark2_palette[1], "RND" = dark2_palette[2]),
    labels = c("DUP" = "True duplicate", "RND" = "False duplicate")
  ) +
  labs(
    x = "Structural divergence ratio",
    y = "Density",
    fill = "Group type"
  ) +
  theme(
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    legend.margin = margin(t = -6)
  )
