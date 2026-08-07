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

inter_gene_cor <- 0.05

mu <- rep_len(0, set_size)

Sigma <- matrix(inter_gene_cor, nrow = set_size, ncol = set_size)
diag(Sigma) <- 1

gene_set <- list(
  "set" = genes[seq_len(set_size)]
)


# ssGSEA2.0 (alpha = 0) ----

message("ssGSEA2.0 (alpha = 0)...")

pvals_ssgsea2_0 <- map_dbl(seq_len(n_datasets), function(i) {
  set.seed(i)

  x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
  y <- rnorm(n_genes - set_size)

  stats <- c(x, y)
  names(stats) <- genes

  ssgsea2(
    stats = stats,
    gene_sets = gene_set,
    alpha = 0,
    seed = i,
    time_only = FALSE
  )[["results"]][["p_value"]]
}, .progress = TRUE)


# ssGSEA2.0 (alpha = 1) ----

message("ssGSEA2.0 (alpha = 1)...")

pvals_ssgsea2_1 <- map_dbl(seq_len(n_datasets), function(i) {
  set.seed(i)

  x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
  y <- rnorm(n_genes - set_size)

  stats <- c(x, y)
  names(stats) <- genes

  ssgsea2(
    stats = stats,
    gene_sets = gene_set,
    alpha = 1,
    seed = i,
    time_only = FALSE
  )[["results"]][["p_value"]]
}, .progress = TRUE)


# HPGSEA (alpha = 0) ----

message("HPGSEA (alpha = 0)...")

pvals_hpgsea_0 <- map_dbl(seq_len(n_datasets), function(i) {
  set.seed(i)

  x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
  y <- rnorm(n_genes - set_size)

  stats <- c(x, y)
  names(stats) <- genes

  hpgsea(
    stats = stats,
    gene_sets = gene_set,
    alpha = 0,
    seed = i
  )[["p_value"]]
}, .progress = TRUE)


# HPGSEA (alpha = 1) ----

message("HPGSEA (alpha = 1)...")

pvals_hpgsea_1 <- map_dbl(seq_len(n_datasets), function(i) {
  set.seed(i)

  x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
  y <- rnorm(n_genes - set_size)

  stats <- c(x, y)
  names(stats) <- genes

  hpgsea(
    stats = stats,
    gene_sets = gene_set,
    alpha = 1,
    seed = i
  )[["p_value"]]
}, .progress = TRUE)


# FGSEA (alpha = 0) ----

message("FGSEA (alpha = 0)...")

pvals_fgsea_0 <- map_dbl(seq_len(n_datasets), function(i) {
  set.seed(i)

  x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
  y <- rnorm(n_genes - set_size)

  stats <- c(x, y)
  names(stats) <- genes

  invisible(
    capture.output(
      pval <- fgsea(
        stats = stats,
        pathways = gene_set,
        gseaParam = 0,
        nproc = 1L
      )[["pval"]]
    )
  )

  return(pval)
}, .progress = TRUE)


# FGSEA (alpha = 1) ----

message("FGSEA (alpha = 1)...")

pvals_fgsea_1 <- map_dbl(seq_len(n_datasets), function(i) {
  set.seed(i)

  x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
  y <- rnorm(n_genes - set_size)

  stats <- c(x, y)
  names(stats) <- genes

  invisible(
    capture.output(
      pval <- fgsea(
        stats = stats,
        pathways = gene_set,
        gseaParam = 1,
        nproc = 1L
      )[["pval"]]
    )
  )

  return(pval)
}, .progress = TRUE)


# t-test ----

message("t-test...")

pvals_t <- map_dbl(seq_len(n_datasets), function(i) {
  set.seed(i)

  x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
  y <- rnorm(n_genes - set_size)

  t.test(x = x, y = y, var.equal = TRUE)[["p.value"]]
}, .progress = TRUE)


# Wilcox test ----

message("Wilcox test...")

pvals_wilcox <- map_dbl(seq_len(n_datasets), function(i) {
  set.seed(i)

  x <- MASS::mvrnorm(n = 1L, mu = mu, Sigma = Sigma)
  y <- rnorm(n_genes - set_size)

  wilcox.test(x = x, y = y)[["p.value"]]
}, .progress = TRUE)


# Stack and save results
res <- list(
  "HPGSEA, 0" = pvals_hpgsea_0,
  "HPGSEA, 1" = pvals_hpgsea_1,
  "ssGSEA2.0, 0" = pvals_ssgsea2_0,
  "ssGSEA2.0, 1" = pvals_ssgsea2_1,
  "FGSEA, 0" = pvals_fgsea_0,
  "FGSEA, 1" = pvals_fgsea_1,
  "Wilcox's rank sum test, -" = pvals_wilcox,
  "t-test, -" = pvals_t
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

saveRDS(
  res,
  file = file.path(
    "simulation",
    "data",
    "pvalues_inter_gene_correlation.rds"
  )
)
