generate_data <- function(nGenes,
                          nSamples,
                          minSetSize,
                          maxSetSize,
                          nSets,
                          nperm,
                          alpha)
{
   on.exit(invisible(gc()))

   # Numeric matrix with genes as rows and samples as columns
   set.seed(0)

   n_digits <- floor(log10(nGenes)) + 1L
   genes <- sprintf(paste0("gene%0", n_digits, "d"), seq_len(nGenes))
   samples <- paste0("sample", seq_len(nSamples))

   # Sample values from a standard normal distribution
   X <- rnorm(n = nGenes * nSamples)
   dim(X) <- c(nGenes, nSamples)
   dimnames(X) <- list(genes, samples)

   # List of gene sets
   size_range <- maxSetSize - minSetSize + 1L
   n_reps <- ceiling(nSets / size_range)
   set_sizes <- rep(minSetSize:maxSetSize, times = n_reps)[seq_len(nSets)]
   set_sizes <- sample(set_sizes) # shuffle sizes

   gene_sets <- lapply(set_sizes, function(size_i) {
      sample(genes, size = size_i)
   })
   names(gene_sets) <- paste0("GeneSet_", seq_along(gene_sets))

   # Results
   out <- list("X" = X,
               "gene_sets" = gene_sets)

   return(out)
}
