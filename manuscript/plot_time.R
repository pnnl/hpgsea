suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggh4x)

  source("manuscript/functions/function-create_dir.R")
})

data_dir <- create_data_dir()
figures_dir <- create_figures_dir()

plot_df <- list(
  "fgsea_1" = "time_fgsea_nproc_1.rds",
  "fgsea_0" = "time_fgsea_nproc_0.rds",
  "hpgsea_1" = "time_hpgsea_nthreads_1.rds",
  "hpgsea_0" = "time_hpgsea_nthreads_0.rds"
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
    nproc = sub(".*_", "", method),
    method = sub("_.*", "", method),
    nproc = paste0(
      ifelse(method == "hpgsea",
             "nthreads = ",
             "nproc = "),
      nproc
    ),
    across(
      .cols = c(nproc, method),
      .fns = ~ factor(.x, levels = unique(.x))
    )
  )

plot_df %>%
  summarise(
    .by = c(method, maxSetSize, nSets, nproc),
    mean_time = mean(elapsed_time)
  ) %>%
  mutate(
    nproc = sub(".* ", "", nproc)
  ) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = mean_time
  ) %>%
  mutate(
    ratio = fgsea / hpgsea
  ) %>%
  arrange(nproc, ratio)

p <- ggplot(plot_df, aes(x = x_pos, y = elapsed_time, color = nSets)) +
  geom_point(
    size = 1.5,
    shape = 21,
    stroke = 0.8
  ) +
  # Cover bottom y-axis grid line
  geom_hline(
    yintercept = 0,
    color = "white"
  ) +
  ggh4x::facet_nested_wrap(
    ~ method + nproc,
    nrow = 1L,
    scales = "free_y"
  ) +
  facetted_pos_scales(
    y = list(
      nproc == "nproc = 1" ~
        scale_y_continuous(
          breaks = seq(0, 64, 16),
          limits = c(0, 64),
          expand = expansion(mult = c(4e-3, 0.05))
        ),
      nproc == "nproc = 0" ~
        scale_y_continuous(
          breaks = seq(0, 16, 4),
          limits = c(0, 16),
          expand = expansion(mult = c(4e-3, 0.05))
        ),
      nproc == "nthreads = 1" ~
        scale_y_continuous(
          breaks = seq(0, 16, 4),
          limits = c(0, 16),
          expand = expansion(mult = c(4e-3, 0.05))
        ),
      nproc == "nthreads = 0" ~
        scale_y_continuous(
          breaks = seq(0, 4, 1),
          limits = c(0, 4),
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
      margin = margin(t = 0, b = 2, unit = "pt")
    ),
    panel.spacing.x = unit(2, "pt"),
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
