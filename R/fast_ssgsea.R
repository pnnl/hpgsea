#' @title High-Performance Variant of Pre-Ranked Gene Set Enrichment Analysis
#'
#' @description A high-Performance variant of pre-ranked Gene Set Enrichment
#'   Analysis (GSEA) that is capable of testing gene sets where each gene has an
#'   expected direction of change (Krug \emph{et al.}, 2019).
#'
#' @param stats a numeric vector with genes/molecules as names. Missing values
#'   are allowed. Zeroes may cause some enrichment scores (ES) to be \code{NA}.
#' @param gene_sets a named list of molecular signatures to test; usually gene
#'   sets. Each element of the list must be a character vector. At least some
#'   elements of \code{gene_sets} must be in \code{names(stats)}.
#' @param alpha numeric (\eqn{\geq 0}); the power to which the absolute values
#'   of the entries of \code{stats} will be raised. Affects the ES calculation.
#'   If \code{alpha=0}, all genes/molecules in each set will contribute equally.
#' @param nperm integer (\eqn{\geq 0}); the number of permutations used to
#'   calculate the normalized enrichment scores (NES) and p-values. Between 0
#'   and 2 billion.
#' @param min_size integer (\eqn{\geq 2}); the minimum set size. To be
#'   considered for testing, sets must have at least \code{min_size} elements
#'   with non-missing values in \code{stats}.
#' @param max_size integer or \code{Inf}; the maximum set size. A value of 500
#'   is recommended.
#' @param sort logical; should the results be sorted by p-value? Default is
#'   \code{TRUE}.
#' @param seed integer or \code{NULL}; if \code{NULL} (default), the normalized
#'   enrichment scores and p-values will vary between runs.
#' @param alternative character; the alternative hypothesis. One of
#'   "\code{two.sided}" (default), "\code{less}", or "\code{greater}". The
#'   latter two will perform one-sided tests.
#'
#' @returns A \code{data.frame} with the following columns:
#'
#' \describe{
#'   \item{set}{character; the gene set that was tested.}
#'
#'   \item{set_size}{integer; number of genes in the set with non-missing values
#'   in \code{stats}.}
#'
#'   \item{ES_u}{numeric; only included if \code{gene_sets} is a directional
#'   database. The enrichment score for the elements of the set that are
#'   expected to be up-regulated.}
#'
#'   \item{ES_d}{numeric; only included if \code{gene_sets} is a directional
#'   database. The enrichment score for the elements of the set that are
#'   expected to be down-regulated.}
#'
#'   \item{ES}{numeric; the enrichment score (ES). The area under the running
#'   sum. This differs from the classic definition of the ES, which is the most
#'   extreme value of the cumulative sum. If \code{gene_sets} is a directional
#'   database, the ES is calculated as \code{ES_u - ES_d}, where more positive
#'   values indicate strong agreement between the true and expected directions
#'   of change and more negative values indicate strong disagreement.}
#'
#'   \item{NES}{numeric; normalized enrichment score (NES). The ratio of the ES
#'   to the absolute mean of the permutation ES with the same sign. If
#'   \code{nperm=0}, all NES will be \code{NA}.}
#'
#'   \item{n_same_sign}{integer; the number of permutation ES with the same sign
#'   as the true ES. At most \code{nperm}. If \code{nperm=0}, all values will be
#'   \code{NA}.}
#'
#'   \item{n_as_extreme}{integer; the number of permutation ES that are at least
#'   as extreme as the true ES. At most \code{n_same_sign}. If \code{nperm=0},
#'   all values will be \code{NA}.}
#'
#'   \item{p_value}{numeric; permutation p-value. Calculated as
#'   \code{(n_as_extreme + 1L) / (n_same_sign + 1L)} if
#'   \code{alternative="two.sided"} (default). If \code{nperm=0}, all values
#'   will be \code{NA}.}
#'
#'   \item{adj_p_value}{numeric; Benjamini and Hochberg FDR adjusted p-value.}
#' }
#'
#' @author Tyler Sagendorf
#'
#' @references Krug, K., Mertins, P., Zhang, B., Hornbeck, P., Raju, R., Ahmad,
#'   R., Szucs, M., Mundt, F., Forestier, D., Jane-Valbuena, J., Keshishian, H.,
#'   Gillette, M. A., Tamayo, P., Mesirov, J. P., Jaffe, J. D., Carr, S. A., &
#'   Mani, D. R. (2019). A Curated Resource for Phosphosite-specific Signature
#'   Analysis. \emph{Molecular & cellular proteomics : MCP, 18}(3), 576–593.
#'   doi:\href{https://doi.org/10.1074/mcp.TIR118.000943}{
#'   10.1074/mcp.TIR118.000943}
#'
#' @export fast_ssgsea

fast_ssgsea <- function(stats,
                        gene_sets,
                        alpha = 1,
                        nperm = 1e5L,
                        min_size = 2L,
                        max_size = Inf,
                        sort = TRUE,
                        seed = NULL,
                        alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(
    arg = alternative,
    choices = c("two.sided", "less", "greater")
  )

  stats <- .prepare_stats(stats)
  n_genes <- length(stats)

  .validate_params(
    alpha = alpha,
    nperm = nperm,
    min_size = min_size,
    max_size = max_size,
    sort = sort,
    seed = seed,
    n_genes = n_genes
  )

  ES_list <- .calc_ES(
    stats = stats,
    alpha = alpha,
    n_genes = n_genes,
    gene_sets = gene_sets,
    min_size = min_size,
    max_size = max_size
  )

  tab <- .calc_ES_perm(
    seed = seed,
    nperm = nperm,
    y = y,
    r = r,
    n_genes = n_genes,
    ES_list = ES_list
  )

  tab <- .calc_pvals(
    tab = tab,
    nperm = nperm,
    sort = sort,
    alternative = alternative
  )

  return(tab)
}
