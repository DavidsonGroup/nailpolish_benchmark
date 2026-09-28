# fmt: skip file

# individual steps of the pipeline and their runtimes
perf_steps <- data.frame(
  method_label = c("nailpolish", "nailpolish", "sicelore", "umitools"),
  step_name = c("np_index", "np_consensus_slim", "sicelore_compute_consensus", "umi_tools_bedtools_dedup")
)

perf_method_labels <- c(nailpolish = "Nailpolish", sicelore = "SiCeLoRe", umitools_bedtools = "UMI-tools")

# ========================
# preprocess data
# ========================

dup_counts <- df.counts %>%
  select(sample, dup_count = pre_dup_groups)

df.perf <- runtime.raw %>%
  inner_join(perf_steps, by = c("process" = "step_name")) %>%
  group_by(sample, method_label) %>%
  summarise(
    runtime_min = sum(s) / 60,
    max_rss_gb = max(max_rss) / 1000,
    .groups = "drop"
  ) %>%
  mutate(
    method_label = factor(method_label, labels = perf_method_labels)
  ) %>%
  left_join(dup_counts, by = "sample")

# write to csv
df.perf %>%
  pivot_wider(
    names_from = method_label,
    values_from = c(runtime_min, max_rss_gb),
    names_sep = "."
  ) %>%
  write_csv("results/tables/performance.csv")


# ========================
# create plot
# ========================

fig.runtime <- df.perf %>%
  ggplot(aes(
    x = dup_count,
    y = runtime_min,
    colour = method_label,
    shape = method_label
  )) +
  geom_point(size = 2) +
  scale_shape_manual(values = method_shapes) +
  scale_colour_manual(values = method_colours) +
  scale_x_log10(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_y_log10() +
  labs(
    x = "Number of input duplicate groups",
    y = "Runtime\n(mins)",
    colour = NULL,
    shape = NULL
  ) +
  theme(
    plot.tag.position = c(-0.18, 1.02),
    plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )


fig.peak_mem <- df.perf %>%
  ggplot(aes(
    x = dup_count,
    y = max_rss_gb,
    colour = method_label,
    shape = method_label
  )) +
  geom_point(size = 2) +
  scale_shape_manual(values = method_shapes) +
  scale_colour_manual(values = method_colours) +
  scale_x_log10(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_y_log10() +
  labs(
    x = "Number of input duplicate groups",
    y = "Peak memory\n(GB)",
    colour = NULL,
    shape = NULL
  ) +
  theme(
    plot.tag.position = c(-0.18, 1.07),
    plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )
