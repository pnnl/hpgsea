library(ggplot2)

source("manuscript/functions/function-create_dir.R")

data_dir <- create_data_dir()
figures_dir <- create_figures_dir()

df <- readRDS(file.path(data_dir, "power.rds"))

levels(df$set_size) <- paste(levels(df$set_size), "genes")
levels(df$method) <- sub("_", ", \u03b1 = ", levels(df$method))

p <- ggplot(df, aes(x = prop_DE, y = power, color = method, linetype = mu)) +
  geom_line(
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ set_size,
    nrow = 1L
  ) +
  scale_x_continuous(
    name = "Proportion of differentially expressed genes",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(
      mult = 5e-3
    )
  ) +
  scale_y_continuous(
    name = "Power",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(
      mult = c(5e-3, 0.02)
    )
  ) +
  scale_color_manual(
    name = "Method:",
    values = c("grey60", "black", "#fdbf6f", "#ff7f00"),
    breaks = levels(df$method)
  ) +
  scale_linetype_manual(
    name = "Effect size:",
    values = c("82", "11", "3132", "solid")
  ) +
  # coord_equal() +
  guides(
    color = guide_legend(
      order = 1L,
      nrow = 1L,
      theme = theme(
        legend.key.width = unit(24, "pt")
      )
    ),
    linetype = guide_legend(
      nrow = 1L,
      theme = theme(
        legend.key.width = unit(24, "pt")
      )
    )
  ) +
  theme_bw(
    base_size = 10,
    base_family = "serif"
  ) +
  theme(
    axis.title.x = element_text(
      color = "black",
      size = 10,
      margin = margin(t = 4, b = 0, unit = "pt")
    ),
    axis.title.y = element_text(
      color = "black",
      size = 10
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.text = element_text(
      color = "black",
      size = 8
    ),
    panel.border = element_rect(
      linewidth = 0.5,
      fill = NA,
      color = "black"
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing.x = unit(12, units = "pt"),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 10,
      margin = margin(t = 1, b = 2, unit = "pt")
    ),
    plot.margin = unit(c(0, 6, 0, 1), "pt"),
    legend.key.height = unit(8, "pt"),
    legend.key.spacing.x = unit(10, "pt"),
    legend.margin = margin(t = -4, b = 0, unit = "pt"),
    legend.title = element_text(
      color = "black",
      size = 10,
      hjust = 0
    ),
    legend.text = element_text(
      color = "black",
      size = 8
    ),
    legend.title.position = "left",
    legend.spacing.x = unit(80, "pt"),
    legend.justification = "left",
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.box.margin = margin(l = 10, unit = "pt"),
    legend.location = "plot"
  )

ggsave(
  filename = file.path(figures_dir, "power.pdf"),
  plot = p,
  height = 2.5,
  width = 5,
  device = cairo_pdf # correctly render unicode
)
