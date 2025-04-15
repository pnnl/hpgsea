#' @title Fast Single-Sample Gene Set Enrichment Analysis (ssGSEA)
#'
#' @description Highly optimized version of the ssGSEA algorithm (Barbie
#'   \emph{et al.}, 2009; Krug \emph{et al.}, 2019).
#'
#' @param X a numeric matrix with genes/molecules as row names and one or more
#'   samples, coefficients, or contrasts as column names. Missing values are
#'   allowed.
#' @param gene_sets a named list of molecular signatures to test. Each element
#'   of the list must be a character vector. Most commonly, molecular signatures
#'   will be in the form of gene sets, though it is not a requirement. At least
#'   some elements of \code{gene_sets} must be in \code{rownames(X)}.
#' @param alpha numeric (\eqn{\geq 0}); the power to which the absolute values
#'   of the entries of \code{X} will be raised. Affects the ES calculation. If
#'   \code{alpha=0}, computation time may be significantly reduced, though all
#'   genes/molecules in each set will contribute equally.
#' @param nperm integer (\eqn{\geq 0}); the number of permutations used to
#'   calculate the normalized enrichment scores (NES) and p-values.
#' @param batch_size integer (\eqn{\gg 1} and \eqn{\leq} \code{nperm}); the
#'   maximum number of permutations run as a single batch. Smaller values will
#'   significantly decrease memory usage and slightly increase runtime.
#' @param adjust_globally logical; whether p-values from different columns of
#'   \code{X} should be adjusted together. Default is \code{FALSE}.
#' @param min_size integer; the minimum set size. To be considered for testing,
#'   sets must have at least \code{min_size} elements with non-missing values in
#'   \emph{all} columns of \code{X}. The default value of 2 is the minimum
#'   possible set size required for testing, though higher values tend to
#'   produce more robust results.
#' @param sort logical; should the results for each column of \code{X} be sorted
#'   by p-value? Default is \code{TRUE}.
#' @param seed integer or \code{NULL}; passed to \code{\link[base]{set.seed}}.
#'   If \code{NULL} (default), the normalized enrichment scores and p-values
#'   will vary between runs.
#'
#' @details Unlike the original ssGSEA implementation, p-values are computed as
#'   described by Phipson and Smyth (2010).
#'
#' @returns A \code{data.frame} with the following columns:
#'
#' \describe{
#'   \item{sample}{factor; one of \code{colnames(X)}.}
#'   \item{set}{character; the gene set being tested.}
#'   \item{set_size}{integer; number of genes in the set with non-missing values
#'   in the \code{mat} matrix for a given column.}
#'   \item{ES_u}{numeric; only included if \code{gene_sets} is a directional
#'   database. The enrichment score for the elements that are expected to be
#'   up-regulated.}
#'   \item{ES_d}{numeric; only included if \code{gene_sets} is a directional
#'   database. The enrichment score for the elements that are expected to be
#'   down-regulated.}
#'   \item{ES}{numeric; the enrichment score (ES). The area under the running
#'   sum. If \code{gene_sets} is a directional database, it is calculated as
#'   \code{ES_u - ES_d}, where more positive values indicate better agreement
#'   between the true and expected directions of change and more negative values
#'   indicate worse agreement.}
#'   \item{NES}{numeric; normalized enrichment score (NES). The ratio of the ES
#'   to the absolute mean of the permutation ES with the same sign. If
#'   \code{nperm=0}, all NES will be \code{NA}.}
#'   \item{n_same_sign}{integer; the number of permutation ES with the same sign
#'   as the true ES. At most \code{nperm}. If \code{nperm=0}, all values will be
#'   \code{NA}.}
#'   \item{n_as_extreme}{integer; the number of permutation ES with the same
#'   sign as the true ES that are at least as extreme as the true ES. At most
#'   \code{n_same_sign}. If \code{nperm=0}, all values will be \code{NA}.}
#'   \item{p_value}{numeric; permutation p-value. Calculated as
#'   \code{(n_as_extreme + 1L) / (n_same_sign + 1L)}.}
#'   \item{adj_p_value}{numeric; Benjamini and Hochberg FDR adjusted p-value.}
#' }
#'
#' @author Tyler Sagendorf
#'
#' @references Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S.
#'   E., Dunn, I. F., Schinzel, A. C., Sandy, P., Meylan, E., Scholl, C.,
#'   Fröhling, S., Chan, E. M., Sos, M. L., Michel, K., Mermel, C., Silver, S.
#'   J., Weir, B. A., Reiling, J. H., Sheng, Q., Gupta, P. B., … Hahn, W. C.
#'   (2009). Systematic RNA interference reveals that oncogenic KRAS-driven
#'   cancers require TBK1. \emph{Nature, 462}(7269), 108–112.
#'   doi:\href{https://doi.org/10.1038/nature08460}{10.1038/nature08460}
#'
#'   Phipson, B., and Smyth, G. K. (2010). Permutation \emph{p}-values should
#'   never be zero: calculating exact p-values when permutations are randomly
#'   drawn. \emph{Stat. Appl. Genet. Molec. Biol.} Volume 9, Issue 1, Article
#'   39.
#'
#'   Krug, K., Mertins, P., Zhang, B., Hornbeck, P., Raju, R., Ahmad, R., Szucs,
#'   M., Mundt, F., Forestier, D., Jane-Valbuena, J., Keshishian, H., Gillette,
#'   M. A., Tamayo, P., Mesirov, J. P., Jaffe, J. D., Carr, S. A., & Mani, D. R.
#'   (2019). A Curated Resource for Phosphosite-specific Signature Analysis.
#'   \emph{Molecular & cellular proteomics : MCP, 18}(3), 576–593.
#'   doi:\href{https://doi.org/10.1074/mcp.TIR118.000943}{
#'   10.1074/mcp.TIR118.000943}
#'
#'   Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., &
#'   Sergushichev, A. (2021). Fast gene set enrichment analysis. \emph{bioRxiv},
#'   060012. doi:\href{https://doi.org/10.1101/060012}{10.1101/060012}
#'
#' @export fast_ssgsea

fast_ssgsea <- function(X,
                        gene_sets,
                        alpha = 1,
                        nperm = 1000L,
                        batch_size = 1000L,
                        adjust_globally = FALSE,
                        min_size = 2L,
                        sort = TRUE,
                        seed = NULL)
{
   # Validate X, sort genes alphabetically, and transpose
   X <- .prepareX(X)

   # Validate function parameters
   .validateParams(alpha = alpha,
                   nperm = nperm,
                   batch_size = batch_size,
                   adjust_globally = adjust_globally,
                   min_size = min_size,
                   sort = sort,
                   seed = seed,
                   n_genes = ncol(X))

   # List of one or two sparse incidence matrices. Genes (rows) are sorted
   # alphabetically.
   A_list <- .sparseIncidence(x = gene_sets,
                              background = colnames(X))

   A <- A_list[["A"]]
   A_d <- A_list[["A_d"]]

   Y <- abs(X) ^ alpha
   R <- .calcRankMatrix(X = X)

   # Avoid propagating NA's when multiplying matrices later
   Z <- !is.na(X)
   storage.mode(Z) <- "numeric"

   NA_idx <- which(Z == 0)
   Y[NA_idx] <- 0
   R[NA_idx] <- 0

   n <- rowSums(Z)

   # Sum of the ranks of entries with nonmissing values
   sumRanks <- n * (n + 1L) / 2L # vector of triangular numbers

   # Calculate set size matrices and remove extreme sets
   M_list <- .calcSetSize(n = n,
                          Z_prime = Z[, rownames(A), drop = FALSE], # Z'
                          A = A,
                          A_d = A_d,
                          min_size = min_size)

   # Extract list components: M, W, M_d, W_d, A (optional), A_d (optional)
   for (name_i in names(M_list))
      assign(name_i, value = M_list[[name_i]])

   # Enrichment score matrices with gene sets as rows and samples as columns
   ES_list <- .calcES(alpha = alpha,
                      Y_prime = Y[, rownames(A), drop = FALSE], # Y'
                      R_prime = R[, rownames(A), drop = FALSE], # R'
                      sumRanks = sumRanks,
                      A = A,
                      M = M,
                      W = W,
                      A_d = A_d,
                      M_d = M_d,
                      W_d = W_d,
                      min_size = min_size)

   ES <- ES_list[["ES"]]
   ES_u <- ES_list[["ES_u"]]
   ES_d <- ES_list[["ES_d"]]

   # Permutations are run in batches to avoid initializing a matrix with nperm
   # columns all at once.
   seed_list <- .createSeedList(nperm = nperm,
                                batch_size = batch_size,
                                seed = seed)

   # List of results for each sample
   tab <- lapply(seq_len(nrow(X)), function(i) {
      # Calculate permutation ES and generate table of results
      tab_i <- .makeResultsTable(
         alpha = alpha,
         nperm = nperm,
         min_size = min_size,
         seed_list = seed_list,
         y_i = Y[i, ],
         r_i = R[i, ],
         n_i = n[i],
         sumRanks_i = sumRanks[i],
         m_i = M[i, ],
         m_d_i = M_d[i, ], # may be NULL
         sets = colnames(A),
         ES_i = ES[i, ],
         # These may be NULL
         ES_u_i = ES_u[i, ],
         ES_d_i = ES_d[i, ]
      )

      return(tab_i)
   })

   names(tab) <- rownames(X)
   tab <- .stackResults(tab = tab,
                        nperm = nperm,
                        sort = sort,
                        adjust_globally = adjust_globally)

   return(tab)
}
