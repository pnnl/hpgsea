library(dqrng)

get_set_sizes <- function(min_size, max_size, n_sets) {
  b <- 10

  y <- log(b + min_size) / log(b + (min_size:max_size))
  y <- (y - min(y)) / diff(range(y))
  lower <- 0.05
  y <- (y * (1 - lower) + lower)

  y <- floor(y / sum(y) * n_sets)

  r <- n_sets - sum(y)

  set_sizes <- rep(min_size:max_size, y)

  set_sizes <- c(set_sizes, rep(min_size:max_size, length.out = r))

  return(set_sizes)
}


generate_data <- function(nGenes,
                          minSetSize,
                          maxSetSize,
                          nSets,
                          nperm) {
  on.exit(invisible(gc()))

  n_digits <- floor(log10(nGenes)) + 1L
  genes <- sprintf(paste0("gene%0", n_digits, "d"), seq_len(nGenes))

  # Gene-level values
  set.seed(0)
  stats <- structure(rnorm(n = nGenes), names = genes)

  # List of gene sets
  size_range <- maxSetSize - minSetSize + 1
  set_sizes <- get_set_sizes(minSetSize, maxSetSize, nSets)

  dqset.seed(0L)
  set_sizes <- dqsample(set_sizes)

  gene_sets <- lapply(set_sizes, function(size_i) {
    dqsample(genes, size = size_i)
  })
  names(gene_sets) <- paste0("GeneSet_", seq_along(gene_sets))

  # Results
  out <- list(
    "stats" = stats,
    "gene_sets" = gene_sets
  )

  return(out)
}
