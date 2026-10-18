suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggh4x)

  source("simulation/functions/function-create_dir.R")
})

data_dir <- create_data_dir()
figures_dir <- create_figures_dir()

plot_df <- list(
  "fgsea_1" = "time_fgsea_nproc_1.rds",
  "fgsea_0" = "time_fgsea_nproc_0.rds",
  "hpgsea" = "time_hpgsea.rds"
) %>%
  lapply(function(f) {
    readRDS(
      file.path(data_dir, f)
    )
  }) %>%
  bind_rows(
    .id = "method"
  ) %>%
  mutate(
    x_pos = 0.2 * (replicate - 2L) + as.numeric(as.factor(maxSetSize)),
    nSets = factor(
      x = nSets,
      levels = sort(unique(nSets)),
      labels = scales::label_comma()(sort(unique(nSets)))
    ),
    method = factor(
      x = method,
      levels = unique(method)
    )
  )

custom_labeller <- function(labels) {
  lapply(labels, function(label) {
    if (label == "hpgsea") {
      bquote("hpgsea" * "," ~ "nperm" == 10^6)
    } else if (label == "fgsea_0") {
      bquote("fgsea" * "," ~ "nproc" == 0)
    } else if (label == "fgsea_1") {
      bquote("fgsea" * "," ~ "nproc" == 1)
    } else {
      stop("Invalid label")
    }
  })
}

p <- ggplot(plot_df, aes(x = x_pos, y = elapsed_time, color = nSets)) +
  geom_point(
    size = 1.8,
    shape = 21,
    stroke = 0.8
  ) +
  # Cover bottom y-axis grid line
  geom_hline(
    yintercept = 0,
    color = "white"
  ) +
  facet_wrap(
    ~ method,
    nrow = 1L,
    scales = "free_y",
    labeller = as_labeller(custom_labeller, default = label_parsed)
  ) +
  facetted_pos_scales(
    y = list(
      # fgsea, nproc = 1
      method == levels(method)[1L] ~ scale_y_continuous(
        breaks = seq(0, 64, 16),
        limits = c(0, 64),
        expand = expansion(mult = c(4e-3, 0.05))
      ),
      # fgsea, nproc = 0; hpgsea
      TRUE ~ scale_y_continuous(
        breaks = seq(0, 16, 4),
        limits = c(0, 16),
        expand = expansion(mult = c(4e-3, 0.05))
      )
    )
  ) +
  scale_x_continuous(
    breaks = 1:2,
    labels = c("500", "1,000"),
    expand = expansion(add = 0.15)
  ) +
  scale_color_manual(
    name = "Number of gene sets:",
    values = c("black", "grey", "#fdbf6f")
  ) +
  labs(
    x = "Maximum gene set size",
    y = "Elapsed time (seconds)"
  ) +
  theme_bw(
    base_family = "serif"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      linetype = "dashed"
    ),
    axis.text = element_text(
      color = "black",
      size = 8
    ),
    axis.title = element_text(
      color = "black",
      size = 10
    ),
    axis.title.x = element_text(
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      margin = margin(r = 2, l = 1)
    ),
    legend.text = element_text(
      color = "black",
      size = 8
    ),
    legend.title = element_text(
      color = "black",
      size = 10
    ),
    legend.justification = "left",
    legend.margin = margin(t = -6, b = 0, l = 10, unit = "pt"),
    legend.location = "plot",
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      color = "black",
      size = 10,
      margin = margin(t = 1, b = 2, unit = "pt")
    ),
    panel.spacing.x = unit(8, "pt"),
    plot.margin = unit(c(0, 1, -0, 0), "pt")
  )

ggsave(
  file.path(figures_dir, "plot_time.pdf"),
  plot = p,
  height = 2.5,
  width = 5,
  bg = "white",
  dpi = 300
)
