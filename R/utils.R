#' @title Validate and Prepare Matrix X
#'
#' @description Validate matrix \code{X}, transpose it, and sort genes
#'   alphabetically. Each column must have at least 3 nonmissing values.
#'
#' @inheritParams fast_ssgsea
#'
#' @returns The transpose of \code{X} with genes sorted alphabetically.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.prepareX <- function(X)
{
   if (!is.matrix(X) ||
       !storage.mode(X) %in% c("integer", "double") ||
       is.null(colnames(X)) ||
       is.null(colnames(X)))
      stop("`X` must be a numeric matrix with row and column names.")

   if (nrow(X) < 3L)
      stop("Matrix `X` must have at least 3 rows.")

   if (any(colSums(!is.na(X)) < 3L))
      stop("Matrix `X` must have at least 3 nonmissing values in each column.")

   # Sort genes alphabetically to deal with ties later
   X <- X[sort(rownames(X)), , drop = FALSE]
   storage.mode(X) <- "numeric" # necessary for certain C++ functions

   X <- t(X)

   return(X)
}


#' @title Validate fast_ssgsea function parameters
#'
#' @inheritParams fast_ssgsea
#'
#' @returns Nothing.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.validateParams <- function(alpha = 1,
                            nperm = 1000L,
                            batch_size = 1000L,
                            adjust_globally = FALSE,
                            min_size = 2L,
                            sort = TRUE,
                            seed = NULL,
                            n_genes)
{
   if (!is.vector(alpha, mode = "numeric") ||
       length(alpha) != 1L ||
       alpha < 0 ||
       is.na(alpha) ||
       is.infinite(alpha))
      stop("`alpha` must be a single non-negative real number.")

   if (!is.vector(nperm, mode = "numeric") ||
       length(nperm) != 1L ||
       is.na(nperm) ||
       nperm < 0L ||
       nperm > 1e6L || # arbitrary limit on number of permutations
       nperm %% 1 != 0) # decimal number
      stop("`nperm` must be a whole number between 0 and 1 million.")

   # batch_size gets modified later, outside of this function
   if (!is.vector(batch_size, mode = "numeric") ||
       length(batch_size) != 1L ||
       is.na(batch_size) ||
       batch_size < min(nperm, 1) ||
       batch_size %% 1 != 0)
      stop("`batch_size` must be a whole number between 1 and `nperm`.")

   batch_size <- min(batch_size, nperm)

   set.seed(seed) # let set.seed validate the seed
   set.seed(NULL)

   if (!is.vector(min_size, mode = "numeric") ||
       length(min_size) != 1L ||
       is.na(min_size) ||
       min_size < 2L ||
       min_size >= n_genes ||
       min_size %% 1 != 0)
      stop("`min_size` must be >= 2 and < nrow(X).")

   if (!is.vector(adjust_globally, mode = "logical") ||
       length(adjust_globally) != 1L ||
       is.na(adjust_globally))
      stop("`adjust_globally` must be TRUE or FALSE.")

   if (!is.vector(sort, mode = "logical") ||
       length(sort) != 1L ||
       is.na(sort))
      stop("`sort` must be TRUE or FALSE.")
}


#' @title Create Sparse Incidence Matrices
#'
#' @description Create a list of sparse incidence matrices, where the unique
#'   sets are rows and the unique elements are columns.
#'
#' @param gene_sets a named list of sets. Each element of the list must be a
#'   character vector. If any elements have the suffix ";d", they will be
#'   separated from the rest and placed in the "A_d" incidence matrix.
#' @param background character vector of elements used to filter the elements of
#'   \code{gene_sets}.
#'
#' @returns A named list of incidence matrices, each of class \code{"dgCMatrix"}
#'   with unique elements as rows and unique sets as columns. \code{"A"} is the
#'   incidence matrix for all elements and \code{"A_d"} is \code{NULL}. If any
#'   elements ended with ";d", \code{"A"} is the incidence matrix of elements
#'   expected to be "up", and \code{"A_d"} is the incidence matrix of elements
#'   expected to be "down". Both matrices will have the same dimensions and
#'   dimension names.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom data.table data.table := chmatch
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
.sparseIncidence <- function(gene_sets, background)
{
   if (!is.list(gene_sets) || is.null(names(gene_sets)))
      stop("`gene_sets` must be a named list of character vectors.")

   elements <- unlist(gene_sets,
                      recursive = FALSE,
                      use.names = FALSE)

   if (!is.vector(elements, mode = "character"))
      stop("`gene_sets` must be a named list of character vectors.")

   # "elements" is a factor, "sets" is not
   dt <- data.table(elements = elements,
                    stringsAsFactors = TRUE)

   dt[, sets := rep(names(gene_sets),
                    lengths(gene_sets))]

   # Determine which elements are expected to be "down", if any
   dt[, direction_down := grepl(";d$", elements, perl = TRUE)]

   # Strip information about direction of change. This may reduce the number of
   # levels if an element is both "up" and "down": "gene;u" and "gene;d" become
   # "gene".
   levels(dt$elements) <- sub(";[ud]{1}$", "", levels(dt$elements))

   # Do not chain with previous line, since the number of levels may change.
   unique_elements <- levels(dt$elements)

   # Convert to characters to use chmatch()
   dt[, elements := as.character.factor(elements)]

   # Only need to check those elements of the background that overlap with the
   # elements of x. No need to validate background, since it gets validated
   # prior to calling this function.
   unique_elements <- intersect(unique_elements, background)

   # Fast character matching (row indices for sparse matrix)
   dt[, i := chmatch(elements, unique_elements, nomatch = 0L)]

   # Remove elements not in the background
   dt <- subset(dt, subset = i != 0L)

   if (!nrow(dt))
      stop("No elements of `gene_sets` are present in rownames(X).")

   unique_sets <- unique(dt[["sets"]])

   # Column indices for sparse matrix
   dt[, j := chmatch(sets, unique_sets)]

   dim_names <- list(unique_elements, unique_sets)
   dims <- lengths(dim_names)

   # Keep genes expected to be "up"
   dt_u <- subset(dt, subset = direction_down == FALSE)

   # Incidence matrix where a 1 indicates that the element is in the set. If x
   # is a directional database, then A will only contain elements that are
   # expected to be "up".
   A <- sparseMatrix(
      i = dt_u[["i"]],
      j = dt_u[["j"]],
      x = 1,
      dims = dims,
      dimnames = dim_names,
      check = FALSE,
      use.last.ij = FALSE
   )

   # In the unlikely event where an element appears multiple times in the same
   # set, some values of A will be > 1. Replace all values with 1. Could also
   # use the use.last.ij parameter in sparseMatrix(), but this is faster.
   attr(A, which = "x") <- rep(1, length(attr(A, which = "x")))

   A_d <- NULL # default

   if (nrow(dt_u) < nrow(dt)) {
      dt_d <- subset(dt, subset = direction_down == TRUE)

      # Incidence matrix where a 1 indicates that a feature is expected to be
      # down in the set.
      A_d <- sparseMatrix(
         i = dt_d[["i"]],
         j = dt_d[["j"]],
         x = 1,
         dims = dims,
         dimnames = dim_names,
         check = FALSE,
         use.last.ij = FALSE
      )

      attr(A_d, which = "x") <- rep(1, length(attr(A_d, which = "x")))

      # The Hadamard product A * A.d should be a matrix of zeros
      if (length(attr(A * A_d, which = "x")))
         stop("Elements can not be both up (suffix \";u\") and ",
              "down (suffix \";d\") in the same set.")
   }

   out <- list("A" = A,
               "A_d" = A_d)

   return(out)
}


#' @title Calculate the Matrix of Ranks of the Gene-Level Values by Sample
#'
#' @param X numeric matrix with samples as rows and genes as columns.
#'
#' @returns A numeric matrix of ranks of the nonmissing values.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom data.table frank
#'
#' @noRd
.calcRankMatrix <- function(X)
{
   R <- apply(X, 1L, function(sample_i) {
      # Columns of X and rows of A are sorted alphabetically. In the case of
      # ties, the largest rank is assigned to the first gene alphabetically.
      frank(sample_i, ties.method = "last", na.last = "keep")
   }, simplify = FALSE)

   R <- do.call(what = rbind, args = R)
   colnames(R) <- colnames(X) # apply(, 1L, ) drops colnames

   storage.mode(R) <- "numeric"

   return(R)
}


#' @title Calculate Set Sizes and Number of Genes Outside of Each Set
#'
#' @param n vector of the number of genes in each sample with nonmissing values.
#' @param Z_prime binary matrix where each row is a sample and each column is a
#'   gene that is present in at least one of the gene sets being tested. A 1
#'   indicates that the corresponding value of matrix \code{t(X)} is nonmissing,
#'   while a 0 indicates missingness.
#' @param A incidence matrix with genes as rows and sets as columns. Indicates
#'   which genes are in each set, or, if \code{A_d} is not \code{NULL},
#'   indicates which genes are expected to be up-regulated in each set.
#' @param A_d incidence matrix with same dimensions as \code{A} or \code{NULL}.
#'   Indicates which genes are expected to be down-regulated in each sample.
#' @param min_size integer; minimum gene set size required for testing. Default
#'   is 2.
#'
#' @returns A named list with the following components:
#'
#' \describe{
#'   \item{"M"}{sample by gene set matrix containing the number of genes with
#'   non-missing values in each set.}
#'
#'   \item{"W"}{sample by gene set matrix containing the number of genes not in
#'   each set with non-missing values in each sample. Calculated as \code{n -
#'   M}.}
#'
#'   \item{"M_d"}{\code{NULL} or the same as \code{"M"}, but only counts genes
#'   expected to be down-regulated in each sample.}
#'
#'   \item{"W_d"}{\code{NULL} or the same as \code{"W"}, but calculated as
#'   \code{n - M_d}.}
#'
#'   \item{"A"}{incidence matrix \code{A} with extremely small or extremely
#'   large gene sets removed. Not included if no sets were removed.}
#'
#'   \item{"A_d"}{incidence matrix \code{A_d} with extremely small or extremely
#'   large gene sets removed. Not included if no sets were removed.}
#' }
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.calcSetSize <- function(n,
                         Z_prime,
                         A,
                         A_d = NULL,
                         min_size = 2L)
{
   # Matrix with samples as rows and genes as columns. Elements are the number
   # of genes in each set with nonmissing values in the sample.
   M <- .Rcpp_matmult_sparse(Z_prime, A) # Z'A
   M[M < min_size] <- 0L

   # Identify sets with too few or too many genes in at least one sample
   extreme_sets_M <- apply(M, 2L, function(m_j) {
      any(m_j == 0L | m_j == n)
   })

   if (!is.null(A_d)) {
      M_d <- .Rcpp_matmult_sparse(Z_prime, A_d)
      M_d[M_d < min_size] <- 0L

      extreme_sets_M_d <- apply(M_d, 2L, function(m_d_j) {
         any(m_d_j == 0L | m_d_j == n)
      })

      extreme_sets <- which(extreme_sets_M & extreme_sets_M_d)
   } else {
      M_d <- W_d <- NULL

      extreme_sets <- which(extreme_sets_M)
   }

   if (length(extreme_sets)) {
      # If any sets are extreme, check if they are all extreme.
      if (length(extreme_sets) == ncol(A))
         stop("All sets in `gene_sets` contain fewer than `min_size` genes ",
              "with nonmissing values or consist of all genes with nonmissing ",
              "values in at least one sample.")

      # Remove extreme sets
      A <- A[, -extreme_sets, drop = FALSE]
      M <- M[, -extreme_sets, drop = FALSE]

      if (!is.null(A_d)) {
         A_d <- A_d[, -extreme_sets, drop = FALSE]
         M_d <- M_d[, -extreme_sets, drop = FALSE]
      }

      out <- list("A" = A,
                  "A_d" = A_d)
   } else {
      out <- list() # no need to update A and A_d
   }

   # Number of genes not in each set with nonmissing values
   W <- n - M

   if (!is.null(A_d))
      W_d <- n - M_d

   out <- c(list("M" = M,
                 "W" = W,
                 "M_d" = M_d,
                 "W_d" = W_d),
            out)

   return(out)
}


#' @title Calculate Matrix of Enrichment Scores
#'
#' @inheritParams fast_ssgsea
#' @param Y_prime the matrix of absolute values of \code{t(X)} raised to the
#'   power of \code{alpha}.
#' @param R_prime the matrix of row vectors of ranks of each element of
#'   \code{t(X)}. Missing values have been replaced with 0.
#' @param sumRanks vector equal to \code{rowSums(R_prime)}.
#' @param A sparse incidence matrix with genes as rows and gene sets as columns.
#'   Only those entries for genes expected to be up-regulated or lacking an
#'   expected direction of change are 1.
#' @param M matrix containing the number of nonmissing values in each set with
#'   samples as rows and sets as columns.
#' @param W matrix containing the number of nonmissing values not in each set
#'   with the same dimensions as \code{M}.
#' @param A_d,M_d,W_d like \code{A}, \code{M_d}, and \code{W_d}, but only for
#'   those genes expected to be down-regulated. If not testing directional sets,
#'   these will be \code{NULL}.
#'
#' @returns A list containing matrices of enrichment scores with names "ES",
#'   "ES_u", and "ES_d". Each matrix has \code{nrow(A)} gene sets as rows and
#'   \code{ncol(Y_prime)} samples as columns. If the sets are not directional,
#'   the "ES_u" and "ES_d" matrices will be \code{NULL}.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.calcES <- function(alpha = 1,
                    min_size = 2L,
                    Y_prime,
                    R_prime,
                    sumRanks,
                    A,
                    M,
                    W,
                    A_d = NULL,
                    M_d = NULL,
                    W_d = NULL)
{
   # Sample by gene set matrix of enrichment scores
   ES_u <- .Rcpp_calcESCore(alpha,
                            min_size,
                            Y_prime,
                            R_prime,
                            sumRanks,
                            A,
                            M,
                            W)

   if (!is.null(A_d)) { # directional database
      ES_d <- .Rcpp_calcESCore(alpha,
                               min_size,
                               Y_prime,
                               R_prime,
                               sumRanks,
                               A_d,
                               M_d,
                               W_d)

      ES <- ES_u - ES_d

      out <- list("ES" = ES,
                  "ES_u" = ES_u,
                  "ES_d" = ES_d)
   } else {
      out <- list("ES" = ES_u)
   }

   return(out)
}


#' @title Calculate Permutation Enrichment Scores for a Single Sample
#'
#' @inheritParams fast_ssgsea
#' @param element_indices integer vector of the indices of nonmissing values in
#'   the sample.
#' @param seeds_batch integer vector of seeds used to generate the permutations.
#' @param y_i dense single row matrix of absolute gene-level values from the
#'   i-th sample raised to the power of \code{alpha}.
#' @param r_i dense single row matrix of the ranks of the gene-level values for
#'   the i-th sample.
#' @param sumRanks_i sum of the ranks for the i-th sample.
#' @param A_perm dense permutation incidence matrix. The number of rows is the
#'   number of unique gene set sizes for the i-th sample. The number of columns
#'   is \code{length(seeds_batch)}.
#' @param theta_m_i vector of unique number of genes with nonmissing values in
#'   each set for the i-th sample.
#' @param theta_w_i vector of unique number of genes with nonmissing values not
#'   in each set for the i-th sample.
#' @param A_perm_d,theta_m_d_i,theta_w_d_i similar to \code{A_perm},
#'   \code{theta_m_i}, and \code{theta_w_i}, but they only contain information
#'   about genes expected to be "down". These will be \code{NULL} if the gene
#'   set database is not directional.
#'
#' @returns A matrix of permutation ES with \code{nrow(A_perm)} rows and
#'   \code{batch_size_b} columns.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom dqrng dqset.seed dqsample.int
#'
#' @noRd
.calcESPerm <- function(alpha = 1,
                        min_size,
                        element_indices,
                        seeds_batch,
                        y_i,
                        r_i,
                        sumRanks_i,
                        A_perm,
                        theta_m_i,
                        theta_w_i,
                        A_perm_d = NULL,
                        theta_m_d_i = NULL,
                        theta_w_d_i = NULL)
{
   max_set_size <- ncol(A_perm)

   n_elements <- length(element_indices) # number of nonmissing values

   batch_size_b <- length(seeds_batch)

   # Integer matrix of indices ranging from 1 to the number of nonmissing values
   # (n_elements). Each column is a different permutation. The number of rows is
   # equal to the size of the largest gene set.
   perm_indices <- vapply(seq_len(batch_size_b), function(p) {
      dqset.seed(seeds_batch[p])

      dqsample.int(n = n_elements, size = max_set_size)
   }, integer(max_set_size)) # integer matrix

   dim(perm_indices) <- NULL # convert matrix to vector

   # Integer vector of indices of nonmissing values for each permutation
   # (permutations are stacked). Faster than element_indices[perm_indices]
   perm_indices <- .Rcpp_indexIntegerVector(element_indices, perm_indices)

   Y_perm <- .Rcpp_indexNumericVector(y_i, perm_indices)
   R_perm <- .Rcpp_indexIntegerVector(r_i, perm_indices)

   # Convert vectors to matrices. The matrices are populated by column.
   dim(Y_perm) <- dim(R_perm) <- c(max_set_size, batch_size_b)

   ES_perm <- .Rcpp_calcESPermCore(alpha,
                                   Y_perm,
                                   R_perm,
                                   sumRanks_i,
                                   A_perm,
                                   theta_m_i,
                                   theta_w_i)

   if (!is.null(A_perm_d)) { # directional sets
      ES_perm_d <- .Rcpp_calcESPermCore(alpha,
                                        Y_perm,
                                        R_perm,
                                        sumRanks_i,
                                        A_perm_d,
                                        theta_m_d_i,
                                        theta_w_d_i)

      # If the number of genes expected to be "up" or "down" in the set is too
      # small, replace the corresponding permutation ES (the entire row) with
      # zero. The final ES will only depend on the "up" or "down" ES. Gene sets
      # were previously filtered so both halves of the ES will not be set to 0.
      ES_perm[theta_m_i < min_size, ] <- 0
      ES_perm_d[theta_m_d_i < min_size, ] <- 0L

      ES_perm <- ES_perm - ES_perm_d
   }

   return(ES_perm)
}


#' @title Create List of Random Seeds for Each Batch of Permutations
#'
#' @inheritParams fast_ssgsea
#'
#' @returns \code{NULL} or a list of random seeds divided into batches.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.createSeedList <- function(nperm = 1000L,
                            batch_size = 1000L,
                            seed = NULL)
{
   if (nperm != 0L) {
      n_batches <- ceiling(nperm / batch_size)

      batch_sizes <- c(rep(batch_size, n_batches - 1L),
                       nperm - batch_size * (n_batches - 1L))

      batch_id <- rep(seq_len(n_batches), batch_sizes)

      # Seeds for permutations
      set.seed(seed)
      seeds <- sample.int(n = 1e7L, size = nperm)

      seed_list <- split(x = seeds, f = batch_id)
      names(seed_list) <- NULL
   } else {
      seed_list <- NULL
   }

   return(seed_list)
}


#' @title Szudzik Pairing Function for Two Integer Vectors
#'
#' @description A pairing function with 100%% packing efficiency that maps a
#'   pair of non-negative integers (order matters) to a single unique integer
#'   (Szudzik, 2006).
#'
#' @param x,y vectors of non-negative integers with the same length.
#'
#' @returns A vector of non-negative integers that uniquely identifies each pair
#'   of values (\code{x}, \code{y}).
#'
#' @author Tyler Sagendorf
#'
#' @references Szudzik, M. (2006). An Elegant Pairing Function. Wolfram Science
#'   Conference. \url{http://szudzik.com/ElegantPairing.pdf}
#'
#' @noRd
.szudzikPairing <- function(x, y)
{
   # Do not need to validate x and y, since this function is only used with
   # integer vectors by other internal functions. Note that x * x prevents
   # integer overflow, unlike x ^ 2.
   x +  ifelse(x < y,
               y * y,      # y ^ 2 + x
               x * x + y)  # x ^ 2 + x + y
}


#' @title Construct Permutation Incidence Matrices and Other Objects
#'
#' @description Construct a list containing permutation incidence matrices and
#'   other objects.
#'
#' @param n_i the total number of nonmissing values in the i-th sample.
#' @param m_i integer vector containing the number of genes in each set or the
#'   number of genes expected to be up-regulated in each set in the i-th sample.
#' @param m_d_i integer vector or \code{NULL}. The number of genes expected to
#'   be down-regulated in each set in the i-th sample. This vector has the same
#'   length as \code{m_i}.
#'
#' @returns A named list with the following components:
#'
#' \describe{
#'   \item{"rep_idx"}{a vector with length \eqn{\geq} \code{ncol(A_perm)} that
#'   maps each row of \code{A_perm} to the corresponding entry of \code{m_i}.
#'   This is used by \code{.extractPermInfo}.}
#'
#'   \item{"A_perm"}{dense incidence matrix where the number of rows is the
#'   number of unique gene set sizes and the number of columns is the size of
#'   the largest gene set. Indicates which genes to use to calculate the
#'   permutation ES or, if \code{"A_perm_d"} is not \code{NULL}, indicates which
#'   genes are expected to be up-regulated to calculate the permutation ES.}
#'
#'   \item{"theta_m_i"}{vector of unique set sizes. If testing a directional
#'   database, the sizes are not necessarily unique, since two gene sets of the
#'   same total size can have different numbers of up- and down-regulated genes,
#'   so those will be treated as separate entries.}
#'
#'   \item{"theta_w_i"}{vector of the same length as \code{m.i}. Calculated as
#'   \code{n_i - theta_m_i}.}
#'
#'   \item{"A_perm_d"}{dense incidence matrix with the same dimensions as
#'   \code{"A_perm"} or \code{NULL}.}
#'
#'   \item{"theta_m_d_i"}{vector with the same length as \code{m_i} containing
#'   the unique numbers of genes expected to be down-regulated in each set or
#'   \code{NULL}. The sizes are not necessarily unique, since two gene sets of
#'   the same total size can have different numbers of up- and down-regulated
#'   genes, so these will be treated as separate entries.}
#'
#'   \item{"theta_w_d_i"}{vector of the number of genes that are not expected to
#'   be "down" in the set. Calculated as \code{n_i - theta_m_d_i}.} }
#'
#' @details Creation of permutation incidence matrices is based on ideas from
#'   Korotkevich \emph{et al.} (2021).
#'
#' @author Tyler Sagendorf
#'
#' @references Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M.
#'   N., & Sergushichev, A. (2021). Fast gene set enrichment analysis.
#'   \emph{bioRxiv}, 060012. doi:\href{
#'   https://doi.org/10.1101/060012}{10.1101/060012}
#'
#' @noRd
.permIncidenceMatrix <- function(n_i,
                                 m_i,
                                 m_d_i = NULL)
{
   # Permutations only need to be generated for each unique set size.
   if (!is.null(m_d_i)) {
      # Include up and down elements or A_perm and A_perm_d may have
      # different dimensions.
      unique_size_mat <- unique(cbind(m_i, m_d_i))

      # Number of genes expected to be up-regulated in each set. May not be
      # unique.
      theta_m_i <- unique_size_mat[, 1L]

      # Number of genes not expected to be up-regulated. May not be unique.
      theta_w_i <- n_i - theta_m_i

      # Number of genes expected to be down-regulated in each set. May not be
      # unique.
      theta_m_d_i <- unique_size_mat[, 2L]

      # Number of genes not expected to be down-regulated. May not be unique.
      theta_w_d_i <- n_i - theta_m_d_i

      # Combine up and down set sizes to get the total number of genes in each
      # set. May not be unique.
      unique_set_sizes <- theta_m_i + theta_m_d_i

      max_set_size <- max(unique_set_sizes)

      A_perm_d <- .Rcpp_calcAPerm(end = theta_m_d_i,
                                  MAX_SET_SIZE = max_set_size,
                                  check = TRUE)

      A_perm <- .Rcpp_calcAPerm(end = unique_set_sizes,
                                MAX_SET_SIZE = max_set_size,
                                check = FALSE)

      A_perm <- A_perm - A_perm_d

      # Each unique pair of entries in m_i and m_d_i are converted to a unique
      # integer. The same is done for each pair of entries in theta_m_i and
      # theta_m_d_i. Then, a vector is generated that maps each (m_i, m_d_i)
      # pair to the unique (theta_m_i, theta_m_d_i) pairs.
      rep_idx <- match(
         .szudzikPairing(m_i, m_d_i),
         .szudzikPairing(theta_m_i, theta_m_d_i)
      )
   } else {
      # Unique number of genes in each set
      theta_m_i <- unique(m_i)

      # Unique number of genes not in each set
      theta_w_i <- n_i - theta_m_i

      theta_m_d_i <- theta_w_d_i <- NULL

      max_set_size <- max(theta_m_i)

      A_perm_d <- NULL

      A_perm <- .Rcpp_calcAPerm(end = theta_m_i,
                                MAX_SET_SIZE = max_set_size,
                                check = FALSE)

      rep_idx <- match(m_i, theta_m_i)
   }

   out <- list(
      "rep_idx" = rep_idx,
      "A_perm" = A_perm,
      "theta_m_i" = theta_m_i,
      "theta_w_i" = theta_w_i,
      # These may be NULL
      "A_perm_d" = A_perm_d,
      "theta_m_d_i" = theta_m_d_i,
      "theta_w_d_i" = theta_w_d_i
   )

   return(out)
}


#' @title Extract Information from a Permutation Enrichment Score Matrix
#'
#' @description Extract information from a matrix of permutation enrichment
#'   scores run as a single batch.
#'
#' @param ES_ls list of enrichment scores grouped by gene set size.
#' @param ES_perm integer matrix of permutation ES. The number of rows is equal
#'   to the length of \code{ES}, while the number of columns is at most the
#'   total number of permutations: more likely, it is a fraction of the total
#'   number of permutations. See the \code{batch_size} parameter of
#'   \code{\link{fast_ssgsea}} for more details.
#'
#' @returns A \code{data.table} with 3 columns:
#'
#' \describe{
#'   \item{"n_same_sign_b"}{integer; the number of permutation ES in each
#'   row of \code{ES_perm} with the same sign as the corresponding ES in
#'   \code{ES}.}
#'   \item{"n_as_extreme_b"}{integer; the number of permutation ES in
#'   each row of \code{ES_perm} that were at least as extreme as the
#'   corresponding ES in \code{ES}. At most \code{"n_same_sign_b"}.}
#'   \item{"sum_ES_perm_b"}{integer; the sum of the absolute values of the
#'   permutation ES that have the same sign as the corresponding ES in
#'   \code{ES}.}
#' }
#'
#' @author Tyler Sagendorf
#'
#' @importFrom data.table data.table := setorderv rbindlist
#'
#' @noRd
.extractPermInfo <- function(ES_ls,
                             ES_perm)
{
   out <- lapply(seq_along(ES_ls), function(i) {
      ES_i <- ES_ls[[i]]

      ES_perm_i <- ES_perm[i, , drop = TRUE]

      out_i <- .Rcpp_extractPermInfo(ES_i, ES_perm_i) # returns list
      class(out_i) <- "data.table"

      return(out_i)
   })

   out <- rbindlist(out)

   return(out)
}


#' @title Generate ssGSEA Results Table for a Single Sample
#'
#' @inheritParams fast_ssgsea
#' @param seed_list list of random seeds for each batch. Ensures that the
#'   permutation enrichment scores will be reproducible.
#' @param y_i numeric vector of absolute gene-level values raised to the power
#'   of \code{alpha} from the i-th sample.
#' @param r_i numeric vector of ranks of the gene-level values from the i-th
#'   sample.
#' @param n_i integer; number of genes in the i-th sample with nonmissing
#'   values.
#' @param sumRanks_i integer; the sum of \code{r_i}.
#' @param m_i the number of genes with nonmissing values in each set from the
#'   i-th sample.
#' @param m_d_i the number of genes with nonmissing values not in each set from
#'   the i-th sample.
#' @param sets character vector of gene set labels.
#' @param ES_i numeric vector of enrichment scores from the i-th sample.
#' @param ES_u_i numeric vector of enrichment scores for the genes in each set
#'   that are expected to be up-regulated in the i-th sample. If not testing a
#'   directional database, this will be \code{NULL}.
#' @param ES_d_i numeric vector of enrichment scores for the genes in each set
#'   that are expected to be down-regulated in the i-th sample. If not testing a
#'   directional database, this will be \code{NULL}.
#'
#' @returns A \code{data.table} with the following columns:
#'
#' \describe{
#'   \item{set}{character; the gene set being tested.}
#'
#'   \item{set_size}{integer; number of genes in the set with non-missing values
#'   in the \code{X} matrix for a given sample.}
#'
#'   \item{ES_u}{numeric; only included if testing a directional database. The
#'   enrichment score for the elements that are expected to be up-regulated.}
#'
#'   \item{ES_d}{numeric; only included if testing a directional database. The
#'   enrichment score for the elements that are expected to be down-regulated.}
#'
#'   \item{ES}{numeric; the enrichment score (ES). The area under the running
#'   sum. If testing a directional database, it is calculated as
#'   \code{ES_u - ES_d}, where more positive values indicate better agreement
#'   between the true and expected directions of change and more negative values
#'   indicate worse agreement.}
#'
#'   \item{NES}{numeric; normalized enrichment score (NES). The ratio of the ES
#'   to the absolute mean of the permutation ES with the same sign. If
#'   \code{nperm=0}, all NES will be \code{NA}.}
#'
#'   \item{n_same_sign}{integer; the number of permutation ES with the same sign
#'   as the true ES. At most \code{nperm}. If \code{nperm=0}, all values will be
#'   \code{NA}.}
#'
#'   \item{n_as_extreme}{integer; the number of permutation ES with the same
#'   sign as the true ES that are at least as extreme as the true ES. At most
#'   \code{n_same_sign}. If \code{nperm=0}, all values will be \code{NA}.}
#'
#'   \item{p_value}{numeric; permutation P-value. Calculated as
#'   \code{(n_as_extreme + 1L) / (n_same_sign + 1L)}.}
#' }
#'
#' @author Tyler Sagendorf
#'
#' @importFrom data.table data.table := setorderv
#'
#' @noRd
.makeResultsTable <- function(alpha = 1,
                              nperm = 1000L,
                              min_size = 2L,
                              seed_list,
                              y_i,
                              r_i,
                              n_i,
                              sumRanks_i,
                              m_i,
                              m_d_i,
                              sets,
                              ES_i,
                              ES_u_i,
                              ES_d_i)
{
   # This will store the results for a single sample
   tab_i <- data.table(
      set = sets,
      set_size = m_i,
      ES = ES_i,
      # Initialize vectors of 0's. These 3 vectors will be updated using the
      # results from each batch of permutations.
      n_same_sign = rep(0L, length(ES_i)),
      n_as_extreme = rep(0L, length(ES_i)),
      sum_ES_perm = rep(0, length(ES_i)),
      row_order = seq_along(ES_i),
      stringsAsFactors = FALSE
   )

   if (nperm != 0L) {
      ## Incidence matrices and other information for permutations
      A_list_perm <- .permIncidenceMatrix(n_i = n_i,
                                          m_i = m_i,
                                          m_d_i = m_d_i)

      # Extract list components: rep_idx, A_perm, theta_m_i, theta_w_i,
      # A_perm_d, theta_m_d_i, theta_w_d_i
      for (name_i in names(A_list_perm))
         assign(x = name_i, value = A_list_perm[[name_i]])

      tab_i[, rep_idx := rep_idx]

      setorderv(tab_i, cols = c("rep_idx", "ES"), order = c(1L, 1L))

      ES_ls <- split(x = tab_i[["ES"]],
                     f = tab_i[["rep_idx"]])
      names(ES_ls) <- NULL

      # Indices of non-missing values for a particular column
      element_indices <- which(r_i != 0L)

      # Optionally split permutations into batches to reduce memory consumption
      for (b in seq_along(seed_list)) {
         ES_perm <- .calcESPerm(alpha = alpha,
                                min_size = min_size,
                                element_indices = element_indices,
                                seeds_batch = seed_list[[b]],
                                y_i = y_i,
                                r_i = r_i,
                                sumRanks_i = sumRanks_i,
                                A_perm = A_perm,
                                theta_m_i = theta_m_i,
                                theta_w_i = theta_w_i,
                                A_perm_d = A_perm_d,
                                theta_m_d_i = theta_m_d_i,
                                theta_w_d_i = theta_w_d_i)

         perm_dt <- .extractPermInfo(ES_ls = ES_ls,
                                     ES_perm = ES_perm)

         # Update summary vectors
         tab_i[, `:=`(
            n_same_sign = n_same_sign + perm_dt[["n_same_sign_b"]],
            n_as_extreme = n_as_extreme + perm_dt[["n_as_extreme_b"]],
            sum_ES_perm = sum_ES_perm + perm_dt[["sum_ES_perm_b"]]
         )]
      } # end permutation batching
   } # end permutations

   setorderv(tab_i, cols = "row_order", order = 1L) # original row order

   if (!is.null(m_d_i))
      tab_i[, `:=`(
         set_size = set_size + m_d_i,
         ES_u = ES_u_i,
         ES_d = ES_d_i
      )]

   return(tab_i)
}


#' @title Stack ssGSEA Results from Multiple Samples
#'
#' @param tab a named list of \code{data.table} objects, each with columns
#'   "set_size", "ES_u" (optional), "ES_d" (optional), "ES", "NES",
#'   "n_same_sign", and "n_as_extreme". Passed to
#'   \code{\link[data.table]{rbindlist}} where the names of the list will be
#'   used to create a column called \code{"sample"}.
#' @param nperm integer; the number of permutations.
#' @param sort logical; whether to sort rows in descending order by p-value.
#' @param adjust_globally logical; whether to adjust all p-values together.
#'
#' @returns A \code{data.frame} with the following columns:
#'
#' \describe{
#'   \item{sample}{factor; one of \code{colnames(X)}.}
#'
#'   \item{set}{character; the gene set being tested.}
#'
#'   \item{set_size}{integer; number of genes in the set with non-missing values
#'   in the \code{X} matrix for a given sample.}
#'
#'   \item{ES_u}{numeric; only included if testing a directional database. The
#'   enrichment score for the elements that are expected to be up-regulated.}
#'
#'   \item{ES_d}{numeric; only included if testing a directional database. The
#'   enrichment score for the elements that are expected to be down-regulated.}
#'
#'   \item{ES}{numeric; the enrichment score (ES). The area under the running
#'   sum. If testing a directional database, it is calculated as
#'   \code{ES_u - ES_d}, where more positive values indicate better agreement
#'   between the true and expected directions of change and more negative values
#'   indicate worse agreement.}
#'
#'   \item{NES}{numeric; normalized enrichment score (NES). The ratio of the ES
#'   to the absolute mean of the permutation ES with the same sign. If
#'   \code{nperm=0}, all NES will be \code{NA}.}
#'
#'   \item{n_same_sign}{integer; the number of permutation ES with the same sign
#'   as the true ES. At most \code{nperm}. If \code{nperm=0}, all values will be
#'   \code{NA}.}
#'
#'   \item{n_as_extreme}{integer; the number of permutation ES with the same
#'   sign as the true ES that are at least as extreme as the true ES. At most
#'   \code{n_same_sign}. If \code{nperm=0}, all values will be \code{NA}.}
#'
#'   \item{p_value}{numeric; permutation P-value. Calculated as
#'   \code{(n_as_extreme + 1L) / (n_same_sign + 1L)}.}
#'
#'   \item{adj_p_value}{numeric; Benjamini and Hochberg FDR adjusted P-value.}
#' }
#'
#' @author Tyler Sagendorf
#'
#' @importFrom data.table rbindlist := setorderv
#' @importFrom stats p.adjust
#'
#' @noRd
.stackResults <- function(tab,
                          nperm = 1000L,
                          sort = TRUE,
                          adjust_globally = FALSE)
{
   sample_names <- names(tab)
   tab <- rbindlist(tab, idcol = "sample")

   tab[, `:=`(
      sample = factor(sample, levels = sample_names),
      set_size = as.integer(set_size),
      NES = ES / (sum_ES_perm / n_same_sign),
      p_value = (n_as_extreme + 1L) / (n_same_sign + 1L)
   )]

   tab[n_same_sign == 0L, # only happens when nperm is very small
       `:=`(NES = NA_real_,
            p_value = NA_real_)]

   if (nperm == 0L) {
      tab[, `:=`(n_same_sign = NA_integer_,
                 n_as_extreme = NA_integer_,
                 p_value = NA_real_)]
   } else if (sort) {
      setorderv(tab, cols = c("sample", "p_value"), order = c(1L, 1L))
   }

   if (adjust_globally) {
      tab[, adj_p_value := p.adjust(p_value, method = "BH")]
   } else {
      tab[, adj_p_value := p.adjust(p_value, method = "BH"),
          by = sample]
   }

   # Reorder/select columns
   keep_cols <- intersect(
      c("sample", "set", "set_size",
        "ES_u", "ES_d", "ES", "NES",
        "n_same_sign", "n_as_extreme",
        "p_value", "adj_p_value"),
      colnames(tab)
   )

   tab <- tab[, keep_cols, with = FALSE]

   # Convert to data.frame
   tab <- as.data.frame(tab)

   return(tab)
}
