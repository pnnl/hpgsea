suppressPackageStartupMessages({
  library(hpgsea)
  library(fgsea)
  library(dplyr)
  library(purrr)
  library(glue)
  library(float)
})

# Set size limits
min_size <- 10L
max_size <- 500L

out_dir <- file.path(
  "simulation",
  "data",
  "test_605_datasets"
)

organisms <- c("human", "mouse")

res <- map(.x = organisms, .f = function(organism) {

  message(
    glue("Testing {organism} datasets...")
  )

  stats_list <- file.path(
    out_dir,
    glue("datasets_{organism}.rds")
  ) %>%
    readRDS() %>%
    lapply(function(li) {
      out <- dbl(li[["t"]])
      names(out) <- li[["Symbol"]]

      return(out)
    })

  dataset_names <- names(stats_list)

  gene_set_file <- ifelse(
    organism == "human",
    "c5.go.v2026.1.Hs.symbols.rds",
    "m5.go.v2026.1.Mm.symbols.rds"
  )

  gene_sets <- file.path(
    out_dir,
    gene_set_file
  ) %>%
    readRDS()

  methods <- c("hpgsea", "fgsea")

  res_method <- map(.x = methods, .f = function(method) {

    res_alpha <- map(.x = c(0, 1), .f = function(alpha) {

      res_data <- map(stats_list, function(stats) {

        if (method == "hpgsea") {

          tic <- Sys.time()

          df <- hpgsea(
            stats = stats,
            gene_sets = gene_sets,
            alpha = alpha,
            nperm = 1e5L,
            min_size = min_size,
            max_size = max_size,
            seed = 0L
          )

          toc <- Sys.time()

          df <- df %>%
            rename(
              size = set_size,
              pval = p_value,
              padj = adj_p_value
            )

        } else if (method == "fgsea") {

          set.seed(0L)

          invisible(
            capture.output({
              tic <- Sys.time()

              df <- fgsea(
                stats = stats,
                pathways = gene_sets,
                minSize = min_size,
                maxSize = max_size,
                gseaParam = alpha,
                nproc = 0L
              )

              toc <- Sys.time()
            })
          )

        }

        time_diff <- as.numeric(difftime(toc, tic, units = "secs"))

        df <- filter(df, !is.na(padj))

        out_i <- data.frame(
          n_sets = nrow(df),
          median_set_size = median(df$size),
          n_signif = sum(df$padj < 0.05),
          min_pval = min(df$pval),
          time_seconds = time_diff
        )

        return(out_i)

      }, .progress = list(
        name = glue("{method}, alpha = {alpha}"),
        clear = FALSE,
        show_after = 0
      )) %>%
        setNames(dataset_names) %>%
        bind_rows(
          .id = "dataset"
        ) %>%
        mutate(
          dataset = factor(
            x = dataset,
            levels = dataset_names
          ),
          across(
            .cols = c(
              contains("signif"),
              starts_with("n_")
            ),
            .fns = as.integer
          )
        )

    }) %>%
      setNames(c(0, 1)) %>%
      bind_rows(
        .id = "alpha"
      )

  }) %>%
    setNames(methods) %>%
    bind_rows(
      .id = "method"
    )

}) %>% # end map organisms
  setNames(organisms) %>%
  bind_rows(
    .id = "organism"
  ) %>%
  mutate(
    across(
      .cols = c(
        organism,
        method,
        alpha
      ),
      .fns = as.factor
    ),
    column_name = paste(method, alpha, sep = "_")
  ) %>%
  arrange(
    organism,
    method,
    alpha,
    dataset
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(
      dataset,
      organism
    ),
    names_from = column_name,
    values_from = c(
      n_sets,
      n_signif,
      median_set_size,
      min_pval,
      time_seconds
    )
  ) %>%
  as.data.frame()

saveRDS(
  object = res,
  file = file.path(
    out_dir,
    "dataset_analysis_results.rds"
  ),
  compress = "xz"
)
