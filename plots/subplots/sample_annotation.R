# Metadata track describing each sample, drawn above the Figure 2 panels.
# fmt: skip file

# ==== DATA ====

# PLACEHOLDER metadata - replace with the real per-sample values.
df.sample_annotation <- tribble(
    ~sample,         ~Platform,        ~Assay,        ~Target,
    "scmixology2",   "ONT\nR9.4.1",     "Single cell", "Non-targeted",
    "GridION",       "ONT\nR9.4.1",     "Single cell", "Non-targeted",
    "GridION Q20",   "ONT\nR10.4",      "Single cell", "Non-targeted",
    "PromethION",    "ONT\nR9.4.1",     "Single cell", "Non-targeted",
    "LR-split-seq",  "PacBio\nSequel II",  "Split-seq",   "Non-targeted",
    "LB SC ONT",     "ONT\nR10.4.1",     "Single cell", "Non-targeted",
    "LB SC PB",      "PacBio\nRevio",  "Single cell", "Non-targeted",
    "LB bulk H146",  "ONT\nR10.4.1",     "Bulk",        "Non-targeted",
    "Visium",        "ONT\nR9.4.1",     "Spatial",     "Non-targeted",
    "DNM1L",         "ONT\nR10.4",     "Bulk",        "Targeted",
    "RAGE-seq",      "ONT\nR9.4.1",     "Single cell", "Targeted",
  ) %>%
  mutate(sample = factor(sample, levels=samples$label)) %>%
  arrange(sample) %>%
  mutate(x = row_number())

row_levels <- c("Platform", "Assay", "Target")
blocks <- df.sample_annotation %>%
  pivot_longer(all_of(row_levels)) %>%
  mutate(name = factor(name, levels = row_levels %>% rev())) %>%
  arrange(name, x) %>%
  group_by(run_id = consecutive_id(value), name, value) %>%
  summarise(xmin = min(x) - 1, xmax = max(x), .groups = "drop") %>%
  mutate(y = cur_group_id(), .by = name)


# ==== PLOT ====

annotation_colours <- setNames(
  c("#AFDDCF", "#F2C7A6", "#F7B4D6"),
  row_levels
)

(fig.sample_annotation <- blocks %>%
  ggplot() +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = y-0.4, ymax = y + 0.4, fill = name),
    colour = "black",
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  geom_text(aes(x = (xmax + xmin) / 2, y = y, label = value), size = 2.2) +
  scale_fill_manual(values = annotation_colours) +
  scale_x_continuous(
    position = "top",
    breaks = df.sample_annotation$x - 0.5,
    labels = df.sample_annotation$sample,
    # expand = c(0.005, 0.005)
    expand = c(0.008, 0.008)
  ) +
  scale_y_continuous(
    breaks = 3:1,
    minor_breaks = NULL,
    labels = row_levels,
    limits = c(0.35, length(row_levels) + 0.65),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x.top = element_text(angle = 15, size = 7, vjust = 0.2),
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#bbbbbb", linewidth = 0.3),
    panel.border = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    plot.margin = margin(t = 4, r = 5.5, b = 0, l = 5.5),
    plot.tag.position = c(-0.08, 1.05)
  ))
