## Given a matrix of samples and vector of genes from a particular set,
## calculate enrichment scores. Used to verify the results of fast_ssgsea.
calculate_ES <- function(X, gene_set, alpha = 1) {
   # Note: only functions imported from data.table for functions in R/ will work
   # without needing `data.table::`
   dt <- data.table::as.data.table(X, keep.rownames = "gene")
   dt[, gene := factor(gene, levels = sort(gene))]

   dt <- data.table::melt(dt,
                          id.vars = "gene",
                          variable.name = "sample",
                          value.name = "x")
   setorderv(dt, cols = c("sample", "x", "gene"), order = c(1L, -1L, 1L))

   dt[, `:=`(y = abs(x) ^ alpha,
             in_set = gene %in% gene_set)]

   # Add column for contribution vector (vector c from the manuscript)
   dt[, c_i := 0]

   # For each gene in the set, divide the y value by the sum of all y values for
   # genes in the set. The result is a value in [0, 1].
   dt[in_set == TRUE, c_i := y / sum(y), by = "sample"]

   # For each gene not in the set, the contribution is negative and only depends
   # on the number of genes not in the set.
   dt[in_set == FALSE, c_i := -1 / .N, by = "sample"]

   # Calculate the cumulative sum of the contributions
   dt[, running_sum := cumsum(c_i), by = "sample"]

   # The sum of the cumulative sum is the ES
   dt[, .(ES = sum(running_sum)), by = "sample"]
}


## Simulate data ----
n_genes <- 2000L
n_samples <- 4L

# Generate values from a standard Normal distribution for testing. Values are
# rounded to the nearest thousandth to introduce ties.
set.seed(0)
X <- matrix(round(rnorm(n_genes * n_samples), digits = 3L),
            nrow = n_genes, ncol = n_samples)

# Sort rows of X by the values in sample 1. This makes it easier to construct
# gene sets with extreme enrichment scores in that sample.
X <- X[order(X[, 1L, drop = TRUE], decreasing = TRUE), , drop = FALSE]

n_digits <- floor(log10(n_genes)) + 1L
gene_names <- sprintf(paste0("gene%0", n_digits, "d"), seq_len(n_genes))

dimnames(X) <- list(gene_names,
                    paste0("sample", seq_len(n_samples)))

# Set1 randomly samples genes that have large positive values in sample 1, while
# Set2 randomly samples genes that have large negative values in sample 2.
gene_set1 <- sample(rownames(X)[1:100], size = 30L)
gene_set2 <- sample(rownames(X)[n_genes - 1:100], size = 60L)

gene_sets <- list("Set1" = gene_set1,
                  "Set2" = gene_set2)
