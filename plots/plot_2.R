# fmt: skip file

source("plots/preprocess.R")

source("plots/subplots/sample_annotation.R")
source("plots/subplots/saturation.R")
source("plots/subplots/median_err_by_sample.R")
source("plots/subplots/performance.R")
source("plots/subplots/fusions.R")


options(repr.plot.width = 6, repr.plot.height = 7)

# The tighter y-axis title margin only applies to the panels that share a
# left-edge alignment; peak memory sits on the right and keeps the default.
tight_y_title <- theme(axis.title.y = element_text(vjust = 0, margin = margin(r = -4)))

# Runtime and memory sit side by side and share one legend.
fig_de <- (
    (fig.runtime / fig.peak_mem) & 
    tight_y_title &
    theme(axis.title.x = element_text(margin = margin(t = -4)))
  ) +
  plot_layout(
    guides = "collect",
    axes = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(0, 0, 3, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(8, "pt")
  )

bottom_panel <- (
  (fig_de | (p_fus_thr + theme(plot.tag.position = c(-0.10, 1.01)))) + 
  plot_layout(widths = c(1, 1))
)

# The annotation track is a flat row rather than a nested plot (patchwork
# flattens nested `/` operands), so it takes an empty tag and the letters carry
# on from the saturation panel.
fig_2 <- (
  ((fig.sample_annotation / fig.saturation / fig.median_err) & tight_y_title) / 
  bottom_panel
) +
  plot_layout(heights = c(0.5, 0.5, 1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F"))) &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.location = "panel", # anchors to panel, not plot margin
    legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(1, 0, 3, 0)
  )

width <- 170
height <- 190

ggsave("results/figures/fig_2.png", fig_2, units = "mm", width = width, height = height, dpi = 300)
ggsave("results/figures/fig_2.pdf", fig_2, units = "mm", width = width, height = height)


sup_fig_2 <- fig.saturation_fdd
ggsave("results/figures/sup_fig_2.png", sup_fig_2, units = "mm", width = 120, height = 80, dpi = 300)
ggsave("results/figures/sup_fig_2.pdf", sup_fig_2, units = "mm", width = 120, height = 80)

sup_fig_4 <- (p_snv_all / p_snv_thr) + plot_layout(
    guides = "collect"
  ) & theme(
    legend.position = "bottom"
  )
ggsave("results/figures/sup_fig_4.png", sup_fig_4, units = "mm", width = 120, height = 135, dpi = 300)
ggsave("results/figures/sup_fig_4.pdf", sup_fig_4, units = "mm", width = 120, height = 135)

options(positron.plot.device = "svg")
