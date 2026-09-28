# fmt: skip file

# ========================
# preprocess data
# ========================

# lookup table with each method and column to access
labels <- data.frame(
  method_label = c("Nailpolish (with FDD)", "Nailpolish (without FDD)", "SiCeLoRe", "UMI-tools", "Random original", "Best original"),
  method = c("nailpolish", "nailpolish_no_cluster", "sicelore", "umitools_bedtools", "nailpolish", "nailpolish"),
  group = c("np_multiread_consensus", "np_multiread_consensus", "duplex", "duplicate", "np_multiread_random_original", "np_multiread_best_original")
) %>%
  mutate(method_label = factor(method_label, levels=rev(method_label)))

df.median_err <- df %>%
  inner_join(labels, by = c("method", "group")) %>%
  arrange(sample, method_label)

# create csv
df.median_err.wide <- df.median_err %>%
  select(sample, method_label, error_rate_median) %>%
  pivot_wider(
    names_from = method_label,
    values_from = error_rate_median
  ) %>%
  mutate(
    `NP vs Best original`   = `Best original`   / `Nailpolish (with FDD)`,
    `NP vs Random original` = `Random original` / `Nailpolish (with FDD)`,
    `NP vs UMI-tools`       = `UMI-tools`       / `Nailpolish (with FDD)`,
    `NP vs SiCeLoRe`        = `SiCeLoRe`        / `Nailpolish (with FDD)`
  )

df.median_err.wide %>%
  write.csv("results/tables/median_err_by_sample.csv")


# ========================
# create plot
# ========================

fig.median_err <- df.median_err %>%
  ggplot(aes(
    x = sample,
    y = error_rate_median,
    colour = method_label,
    fill = method_label,
    shape = method_label
  )) +
  geom_point(size = 3, position = position_dodge(width = 0.6, reverse = TRUE)) +
  scale_shape_manual(values = method_shapes) +
  scale_colour_manual(values = method_colours) +
  scale_fill_manual(values = method_colours) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    x = NULL,
    y = "Median error rate\nof duplicate groups",
    colour = NULL,
    fill = NULL,
    shape = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.27, 0.75),
    plot.tag.position = c(-0.08, 0.95),
    plot.margin = margin(t = 4, r = 5.5, b = 10, l = 5.5)
  ) +
  guides(
    colour = guide_legend(nrow = 3, reverse = TRUE),
    fill = guide_legend(nrow = 3, reverse = TRUE),
    shape = guide_legend(nrow = 3, reverse = TRUE)
  )
