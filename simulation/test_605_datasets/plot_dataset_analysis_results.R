suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(ggplot2)
  library(GGally)
  library(grid)
})

out_dir <- file.path(
  "simulation",
  "figures"
)

df <- file.path(
  "simulation",
  "data",
  "test_605_datasets",
  "dataset_analysis_results.rds"
) %>%
  readRDS()

make_pairs_plot <- function(df) {

  lower_fun <- function(data, mapping) {
    limits <- c(0, 3000)
    breaks <- seq(0, 3000, 1000)

    ggplot(data = data, mapping = mapping) +
      geom_point(
        shape = 16,
        size = 0.7
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        color = "#ff7f00",
        linewidth = 0.5,
        linetype = "53"
      ) +
      scale_x_continuous(
        name = "x-axis label",
        breaks = breaks,
        labels = scales::label_comma(),
        limits = limits
      ) +
      scale_y_continuous(
        breaks = breaks,
        labels = scales::label_comma(),
        limits = limits
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          color = "black",
          size = 8,
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        axis.text.y = element_text(
          color = "black",
          size = 8
        ),
        panel.grid = element_blank(),
        axis.ticks = element_line(
          color = "black",
          linewidth = 0.35
        )
      )
  }

  ggpairs_lower_tri <- function(g) {
    g$plots <- g$plots[-seq_len(g$nrow)]
    g$yAxisLabels <- g$yAxisLabels[-1L]
    g$nrow <- g$nrow - 1L

    g$plots <- g$plots[-seq(g$ncol, length(g$plots), by = g$ncol)]
    g$xAxisLabels <- g$xAxisLabels[-g$ncol]
    g$ncol <- g$ncol - 1L

    return(g)
  }

  prefix <- "n_signif_"

  pair_cols <- paste0(
    prefix,
    c(
      "fgsea_1",
      "fgsea_0",
      "hpgsea_1",
      "hpgsea_0"
    )
  )

  column_labels <- sub(
    "_",
    ", \u03b1 = ", # unicode alpha
    toupper(sub(prefix, "", pair_cols, fixed = TRUE)),
    fixed = TRUE
  )

  axis_title <- "Number of significant gene sets"

  p <- ggpairs(
    data = df,
    columns = pair_cols,
    columnLabels = column_labels,
    lower = list(
      continuous = wrap(lower_fun)
    ),
    upper = "blank",
    diag = list(
      continuous = "blankDiag",
      discrete = "blankDiag"
    ),
    switch = "both",
    showStrips = FALSE,
    axisLabels = "show",
    xlab = axis_title,
    ylab = axis_title
  ) %>%
    ggpairs_lower_tri() +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.35
      ),
      # The PDF won't be centered unless this is done
      plot.margin = unit(c(0, 0, 24, 24), "pt"),
      strip.placement = "outside",
      strip.background = element_rect(
        color = NULL,
        fill = "grey94"
      ),
      strip.switch.pad.grid = unit(18, "pt"),
      strip.text = element_text(
        color = "black",
        size = 10
      ),
      panel.spacing = unit(2, "pt"),
      plot.title = element_text(
        color = "black",
        size = 8,
        margin = margin(t = -0.5, r = 0, b = 0.1, l = -0.5, unit = "in")
      ),
      # The negative margin moves the titles between the facet strips and the
      # axis text
      axis.title.x = element_text(
        color = "black",
        size = 9,
        margin = margin(t = -34, r = 0, b = 0, l = 0, "pt")
      ),
      axis.title.y = element_text(
        color = "black",
        size = 9,
        margin = margin(t = 0, r = -34, b = 0, l = 0, "pt")
      )
    )

  # 1:1 aspect ratio
  p <- ggmatrix_gtable(p)

  panel_size <- unit(1.35, "in")

  p$widths[unitType(p$widths) == "null"] <- panel_size
  p$heights[unitType(p$heights) == "null"] <- panel_size

  return(p)
}

p <- make_pairs_plot(df)

ggsave(
  filename = file.path(
    out_dir,
    "scatterplots_605_datasets.pdf"
  ),
  plot = p,
  height = 5,
  width = 5,
  device = cairo_pdf # render unicode correctly
)
