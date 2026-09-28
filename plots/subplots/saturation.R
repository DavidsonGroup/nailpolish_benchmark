# What fraction of reads in each sample would be removed by consensus calling.
#
# Counts come from data/read_counts.csv. Each duplicate group collapses to a
# single consensus read, so the removed reads are dup_reads - dup_groups. The
# `pre_` columns are before false duplicate detection, the `post_` columns after.

# ==== DATA ====

df.counts %>%
  mutate(saturation_change = post_saturation - pre_saturation) %>%
  write_csv("results/tables/saturation.csv")


# ==== PLOT ====

stage_colours <- c("Without FDD" = dark2_palette[1], "With FDD" = "grey80")

fig.saturation <- df.counts %>%
  ggplot(aes(sample, pre_saturation)) +
  geom_col(width = 0.6, fill = dark2_palette[1]) +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0.02, 0.02)),
    limits = c(0, 1)
  ) +
  labs(x = NULL, y = "Saturation\n(1 - UMIs/reads)", fill = NULL) +
  theme(
    # Sample names live on the annotation track above this panel.
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 4, r = 5.5, b = 4, l = 5.5),
    legend.position = c(0.185, 0.85),
    plot.tag.position = c(-0.08, 1.05)
  ) +
  guides(
    fill = guide_legend(nrow = 1)
  )

# Standalone version comparing saturation before and after false duplicate
# detection, so it carries its own sample labels.
fig.saturation_fdd <- df.counts %>%
  pivot_longer(
    c(pre_saturation, post_saturation),
    names_to = "stage",
    values_to = "saturation"
  ) %>%
  mutate(
    stage = factor(
      stage,
      levels = c("pre_saturation", "post_saturation"),
      labels = c("Without FDD", "With FDD")
    )
  ) %>%
  ggplot(aes(sample, saturation, fill = stage)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = stage_colours) +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0.02, 0.02)),
    limits = c(0, 1)
  ) +
  labs(x = NULL, y = "Saturation\n(1 - UMIs/reads)", fill = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom"
  ) +
  guides(
    fill = guide_legend(nrow = 1)
  )
