suppressPackageStartupMessages({
  library(hpgsea)
  library(fgsea)
  library(dplyr)
  library(purrr)
  library(BiocParallel)
})

n_datasets <- 1e3L
n_genes <- 1e4L
genes <- paste0("gene", seq_len(n_genes))

BPPARAM <- bpparam()

get_pvals <- function(method, set_size, n_DE, mu) {

  gene_set <- list(
    "set" = genes[seq_len(set_size)]
  )

  # Since fgsea also uses bplapply internally, this will constantly trigger the
  # "failed to open port" error. This doesn't cause any issues with the
  # analysis, but it will cause the progress bar to be printed multiple times.
  pvals <- bplapply(
    BPPARAM = BPPARAM,
    X = seq_len(n_datasets),
    FUN = function(i) {
      set.seed(i)

      stats <- rnorm(n_genes)
      stats[seq_len(n_DE)] <- stats[seq_len(n_DE)] + mu
      names(stats) <- genes

      switch(
        method,
        hpgsea_0 = {
          pval <- hpgsea(
            stats = stats,
            gene_sets = gene_set,
            alpha = 0,
            seed = 0L
          )[["p_value"]]
        },
        hpgsea_1 = {
          pval <- hpgsea(
            stats = stats,
            gene_sets = gene_set,
            alpha = 1,
            seed = 0L
          )[["p_value"]]
        },
        fgsea_0 = {
          invisible({
            capture.output({
              set.seed(0L)
              pval <- fgsea(
                pathways = gene_set,
                stats = stats,
                gseaParam = 0,
                nproc = 1L
              )[["pval"]]
            })
          })
        },
        fgsea_1 = {
          invisible({
            capture.output({
              set.seed(0L)
              pval <- fgsea(
                pathways = gene_set,
                stats = stats,
                gseaParam = 1,
                nproc = 1L
              )[["pval"]]
            })
          })
        }
      )

      return(pval)
    }
  ) %>%
    unlist()

  return(pvals)
}

methods <- c("fgsea_0", "fgsea_1", "hpgsea_0", "hpgsea_1")

res <- lapply(methods, function(method) {

  expand.grid(
    set_size = c(20L, 100L, 200L),
    prop_DE = seq(0.1, 1, 0.1),
    mu = c(0.2, 0.5, 1, 2)
  ) %>%
    mutate(
      n_DE = ceiling(prop_DE * set_size),
      pvals = pmap(
        .l = list(method, set_size, n_DE, mu),
        .f = get_pvals,
        .progress = list(
          name = method,
          clear = FALSE,
          show_after = 0
        )
      )
    ) %>%
    tidyr::unnest(
      cols = pvals
    )

}) %>%
  setNames(
    methods
  ) %>%
  bind_rows(
    .id = "method"
  ) %>%
  mutate(
    across(
      .cols = c(
        set_size,
        mu
      ),
      .fns = ~ factor(
        x = .x,
        levels = sort(unique(.x))
      )
    ),
    method = factor(
      x = method,
      levels = methods
    )
  ) %>%
  summarise(
    .by = c(
      method,
      mu,
      prop_DE,
      set_size
    ),
    power = mean(pvals < 0.05)
  )

saveRDS(
  res,
  file = file.path(
    "simulation",
    "data",
    "power_simulation.rds"
  ),
  compress = "xz"
)
