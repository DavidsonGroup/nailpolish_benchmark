# fmt: skip file

library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(fs)


# ============================================================
# CONFIG
# ============================================================

theme_set(
  theme_bw(base_size = 9) +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_line(),
      legend.background = element_blank()
    )
)
options(
  ggplot2.discrete.fill = brewer.pal(8, "Dark2"),
  ggplot2.discrete.colour = brewer.pal(8, "Dark2")
)

# Samples are ordered roughly by error rate, so trends read left to right.
samples <- tribble(
  ~id,                   ~label,
  "lr-split-seq",        "LR-split-seq",
  "longbench_sc_pb",     "LB SC PB",
  "longbench_sc_ont",    "LB SC ONT",
  "gridION_Q20",         "GridION Q20",
  "gridION",             "GridION",
  "scmixology2",         "scmixology2",
  "promethION",          "PromethION",
  "visium",              "Visium",
  "longbench_bulk_H146", "LB bulk H146",
  "dnm1l_barcode01",     "DNM1L",
  "rage-seq",            "RAGE-seq",
)

dark2_palette <- brewer.pal(8, "Dark2")
method_shapes <- c("Nailpolish (with FDD)" = 16, "Nailpolish (without FDD)" = 15, "SiCeLoRe" = 17,
                   "UMI-tools" = 18, "Random original" = 3, "Best original" = 4, "Nailpolish" = 16,
                   "UMI-tools (bedtools)" = 18, "UMI-tools (per cell)" = 9)
method_colours <- c("Nailpolish (with FDD)" = dark2_palette[1], "Nailpolish (without FDD)" = dark2_palette[2],
                    "SiCeLoRe" = dark2_palette[3], "UMI-tools" = dark2_palette[4], "Random original" = dark2_palette[5],
                    "Best original" = dark2_palette[6], "Nailpolish" = dark2_palette[1])

# The performance panels distinguish the two UMI-tools modes, since their cost
# differs substantially; the accuracy panels collapse them to one "UMI-tools".
group_colours <- c("Duplicate" = dark2_palette[1], "Non-duplicate" = "grey80", "No UMI" = dark2_palette[4])

# ============================================================
# Load data
# ============================================================

SAMPLES_TO_DROP = c("longbench_sc_ont_no_chimera", "longbench_bulk_H146_padded", "longbench_bulk_H146_trimmed", "max_2nd")

# Groups not listed here sort after all listed groups, alphabetically.
GROUP_ORDER <- c("np_worst_original", "np_random_original", "np_best_original",
                 "np_consensus", "np_multiread_worst_original", "np_multiread_random_original",
                 "np_multiread_best_original", "np_multiread_consensus", "np_simplex")

# One CSV per method, under data/csv/<sample>/<method>.csv.
df <- dir_ls("data/csv", recurse = TRUE, glob = "*.csv") %>%
  map_dfr(read_csv, show_col_types = FALSE, .id = "path") %>%
  mutate(
    sample = path_file(path_dir(path)),
    method = path_ext_remove(path_file(path)),
    .keep = "unused", .before = 1
  ) %>%
  filter(! sample %in% SAMPLES_TO_DROP) %>%
  mutate(
    sample = factor(sample, levels = samples$id, labels=samples$label),
    group = factor(group, levels = c(GROUP_ORDER, setdiff(sort(unique(group)), GROUP_ORDER))),
    error_rate_mean_pct = error_rate_mean * 100,
    error_rate_median_pct = error_rate_median * 100
  ) %>%
  arrange(sample, method, group)

# load benchmark runtimes
runtime.raw <- dir_ls("data/benchmarks", recurse = TRUE, glob = "*.benchmark.txt") %>%
  map_dfr(read_tsv, show_col_types = FALSE, .id = "path", col_types = "n_nnnnnnnn") %>%
  mutate(
    sample = factor(
      path_file(path_dir(path)),
      levels = samples$id,
      labels = samples$label
    ),
    process = str_remove(path_file(path), "\\.benchmark\\.txt$"),
    .keep = "unused", .before = 1
  )


df.counts <- read_csv("data/read_counts.csv", show_col_types = FALSE) %>%
  filter(! sample %in% SAMPLES_TO_DROP) %>%
  mutate(
    sample = factor(sample, levels = samples$id, labels = samples$label),
    pre_saturation = 1 - (pre_sin + pre_dup_groups) / (pre_sin + pre_dup_reads),
    post_saturation = 1 - (post_sin + post_dup_groups) / (post_sin + post_dup_reads)
  ) %>%
  arrange(sample)

df.counts_old <- df %>%
  filter(method == "nailpolish_no_cluster",
         group %in% c("np_simplex", "np_multiread_consensus")) %>%
  transmute(sample,
            group = factor(group,
                           levels = c("np_simplex", "np_multiread_consensus"),
                           labels = c("singleton", "duplicate")),
            count)
