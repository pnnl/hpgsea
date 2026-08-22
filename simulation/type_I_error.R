suppressPackageStartupMessages({
  library(hpgsea)
  library(fgsea)
  library(dplyr)
  library(purrr)

  source("simulation/scripts/functions/function-ssgsea2.R")
})

n_datasets <- 1e4L # p-values to 4 decimal places
n_genes <- 1e4L
genes <- paste0("gene", seq_len(n_genes))
set_size <- 100L

gene_set <- list(
  "set" = genes[seq_len(set_size)]
)

mu <- rep_len(0, set_size)

for (inter_gene_cor in c(0, 0.05)) {

  message(
    "inter-gene correlation = ", inter_gene_cor, "..."
  )

  Sigma <- matrix(inter_gene_cor, nrow = set_size, ncol = set_size)
  diag(Sigma) <- 1

  get_pvals <- function(method) {

    gene_set <- list(
      "set" = genes[seq_len(set_size)]
    )

    pvals <- map_dbl(.x = seq_len(n_datasets), .f = function(i) {
      set.seed(i)

      if (inter_gene_cor == 0) {
        stats <- rnorm(n_genes)
      } else {
        x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
        y <- rnorm(n_genes - set_size)
        stats <- c(x, y)
      }

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
        ssgsea2_0 = {
          pval <- ssgsea2(
            stats = stats,
            gene_sets = gene_set,
            alpha = 0,
            seed = 0L,
            time_only = FALSE
          )[["results"]][["p_value"]]
        },
        ssgsea2_1 = {
          pval <- ssgsea2(
            stats = stats,
            gene_sets = gene_set,
            alpha = 1,
            seed = 0L,
            time_only = FALSE
          )[["results"]][["p_value"]]
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
    }, .progress = list(
      name = method,
      clear = FALSE,
      show_after = 0
    ))

    return(pvals)
  }

  methods <- c(
    "hpgsea, 0" = "hpgsea_0",
    "hpgsea, 1" = "hpgsea_1",
    "ssGSEA2.0, 0" = "ssgsea2_0",
    "ssGSEA2.0, 1" = "ssgsea2_1",
    "fgsea, 0" = "fgsea_0",
    "fgsea, 1" = "fgsea_1"
  )

  res <- lapply(unname(methods), get_pvals) %>%
    setNames(
      names(methods)
    ) %>%
    tibble::enframe(
      name = "method",
      value = "pval"
    ) %>%
    tidyr::separate_wider_delim(
      cols = method,
      delim = ", ",
      names = c("method", "alpha")
    ) %>%
    mutate(
      across(
        .cols = c(
          method,
          alpha
        ),
        .fns = ~ factor(
          x = .x,
          levels = unique(.x)
        )
      )
    ) %>%
    tidyr::unnest(
      pval
    ) %>%
    summarise(
      .by = c(
        method,
        alpha
      ),
      p0.01 = mean(pval < 0.01),
      p0.02 = mean(pval < 0.02),
      p0.05 = mean(pval < 0.05),
      p0.10 = mean(pval < 0.10)
    )

  if (inter_gene_cor == 0) {
    saveRDS(
      res,
      file = file.path(
        "simulation",
        "data",
        "pvalues_type_I_error.rds"
      )
    )
  } else {
    saveRDS(
      res,
      file = file.path(
        "simulation",
        "data",
        "pvalues_inter_gene_correlation.rds"
      )
    )
  }

}
