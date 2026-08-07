suppressPackageStartupMessages({
  library(hpgsea)
  library(fgsea)
  library(dplyr)
  library(purrr)
})

n_datasets <- 1e3L
n_genes <- 1e4L
genes <- paste0("gene", seq_len(n_genes))

get_pvals <- function(method, set_size, n_DE, mu) {

  method <- match.arg(
    arg = method,
    choices = c("hpgsea", "fgsea", "t-test")
  )

  gene_set <- list(
    "set" = genes[seq_len(set_size)]
  )

  pvals <- map_dbl(.x = seq_len(n_datasets), .f = function(i) {
    set.seed(i)

    stats <- c(
      rnorm(n_DE, mean = mu),
      rnorm(n_genes - n_DE)
    )

    names(stats) <- genes

    switch(
      method,
      hpgsea = {
        pval <- hpgsea(
          stats = stats,
          gene_sets = gene_set,
          seed = 0L
        )[["p_value"]]
      },
      fgsea = {
        invisible({
          capture.output({
            set.seed(0L)
            pval <- fgseaSimple(
              pathways = gene_set,
              stats = stats,
              nperm = 1e3L,
              nproc = 1L,
            )[["pval"]]
          })
        })
      },
      `t-test` = {
        pval <- t.test(
          x = stats[seq_len(set_size)],
          y = stats[(set_size + 1L):n_genes],
          var.equal = TRUE
        )[["p.value"]]
      }
    )

    return(pval)
  })

  return(pvals)
}

methods <- c("hpgsea", "fgsea", "t-test")

res <- lapply(methods, function(method) {
  message(method, "...")

  expand.grid(
    set_size = c(20L, 100L, 200L),
    prop_DE = c(0.1, seq(0.2, 1, 0.2)),
    mu = c(0.25, 0.5, 1, 2)
  ) %>%
    mutate(
      n_DE = ceiling(prop_DE * set_size),
      pvals = pmap(
        .l = list(method, set_size, n_DE, mu),
        .f = get_pvals,
        .progress = TRUE
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
    set_size = factor(
      x = set_size,
      levels = sort(set_size),
      labels = paste0("m = ", sort(set_size))
    ),
    mu = factor(
      x = mu,
      levels = sort(unique(mu))
    ),
    method = factor(
      x = method,
      levels = methods,
      labels = c("HPGSEA", "FGSEA", "t-test")
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
