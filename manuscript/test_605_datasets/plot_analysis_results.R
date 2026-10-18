suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)

  source("manuscript/functions/function-create_dir.R")
})

data_dir <- file.path(create_data_dir(), "test_605_datasets")

if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

figures_dir <- create_figures_dir()

plot_df <- file.path(
  data_dir,
  "analysis_results.rds"
) %>%
  readRDS() %>%
  select(
    dataset,
    starts_with("n_signif")
  ) %>%
  tidyr::pivot_longer(
    cols = contains("hpgsea"),
    names_to = "hpgsea_alpha",
    values_to = "n_signif_hpgsea"
  ) %>%
  tidyr::pivot_longer(
    cols = contains("fgsea"),
    names_to = "fgsea_alpha",
    values_to = "n_signif_fgsea"
  ) %>%
  mutate(
    across(
      .cols = contains("alpha"),
      # u03b1 is unicode alpha
      .fns = ~ sub(".*_([^_]+)_(\\d{1})$", "\\1, \u03b1 = \\2", .x)
    )
  )

axis_title <- "Number of gene sets with P-values < 0.05"
axis_limits <- c(0, 3000)
axis_breaks <- seq(0, 3000, 1000)
point_color <- "#fdbf6f" #"#ff7f00"
text_color <- "grey60"
text_size <- 10 / .pt

n_df <- plot_df %>%
  summarise(
    .by = contains("alpha"),
    # >, not >= to avoid biasing the counts toward hpgsea
    n_above = sum(n_signif_hpgsea > n_signif_fgsea),
    n_below = sum(n_signif_hpgsea < n_signif_fgsea)
  )

p <- ggplot(plot_df, aes(x = n_signif_fgsea, y = n_signif_hpgsea)) +
  geom_point(
    shape = 21,
    size = 1.4,
    stroke = 0.4
  ) +
  geom_text(
    data = n_df,
    aes(x = 500, y = 2500, label = n_above),
    size = text_size,
    color = text_color
  ) +
  geom_text(
    data = n_df,
    aes(x = 2500, y = 500, label = n_below),
    size = text_size,
    color = text_color
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = point_color,
    linewidth = 0.5,
    linetype = "53"
  ) +
  facet_grid(
    rows = vars(hpgsea_alpha),
    cols = vars(fgsea_alpha),
    switch = "both"
  ) +
  scale_x_continuous(
    breaks = axis_breaks,
    labels = scales::label_comma(),
    limits = axis_limits
  ) +
  scale_y_continuous(
    breaks = axis_breaks,
    labels = scales::label_comma(),
    limits = axis_limits
  ) +
  labs(
    x = axis_title,
    y = axis_title
  ) +
  theme_bw(
    base_size = 10,
    base_family = "serif"
  ) +
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    # This is necessary to center the PDF and avoid clipping
    plot.margin = unit(c(1, 1, 19, 19), "pt"),
    strip.placement = "outside",
    strip.background = element_rect(
      color = "grey75",
      fill = "grey94"
    ),
    strip.switch.pad.grid = unit(18, "pt"),
    strip.text = element_text(
      color = "black",
      size = 10,
      margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
    ),
    panel.spacing = unit(4, "pt"),
    # The negative margin moves the titles between the facet strips and the
    # axis text
    axis.title.x = element_text(
      color = "black",
      size = 10,
      margin = margin(t = -28, r = 0, b = 0, l = 0, "pt")
    ),
    axis.title.y = element_text(
      color = "black",
      size = 10,
      margin = margin(t = 0, r = -28, b = 0, l = 0, "pt")
    ),
    axis.text = element_text(
      color = "black",
      size = 8
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.35
    )
  )

ggsave(
  filename = file.path(figures_dir, "scatterplots_605_datasets.pdf"),
  plot = p,
  height = 3.5,
  width = 3.5,
  device = cairo_pdf # render unicode correctly
)
