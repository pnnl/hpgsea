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
.prepareX <- function(X) {
  if (
    !is.matrix(X) ||
      !storage.mode(X) %in% c("integer", "double") ||
      is.null(rownames(X)) ||
      is.null(colnames(X))
  ) {
    stop("`X` must be a numeric matrix with row and column names.")
  }

  if (nrow(X) < 3L) {
    stop("Matrix `X` must have at least 3 rows.")
  }

  if (any(colSums(!is.na(X)) < 3L)) {
    stop("Matrix `X` must have at least 3 nonmissing values in each column.")
  }

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
                            max_size = Inf,
                            sort = TRUE,
                            seed = NULL,
                            n_genes) {
  if (
    !is.vector(alpha, mode = "numeric") ||
      length(alpha) != 1L ||
      alpha < 0 ||
      is.na(alpha) ||
      is.infinite(alpha)
  ) {
    stop("`alpha` must be a single non-negative real number.")
  }

  if (
    !is.vector(nperm, mode = "numeric") ||
      length(nperm) != 1L ||
      is.na(nperm) ||
      nperm < 0L ||
      nperm > 1e6L || # arbitrary limit on number of permutations
      nperm %% 1 != 0 # decimal number
  ) {
    stop("`nperm` must be a whole number between 0 and 1 million.")
  }

  # batch_size gets modified later, outside of this function
  if (
    !is.vector(batch_size, mode = "numeric") ||
      length(batch_size) != 1L ||
      is.na(batch_size) ||
      batch_size < min(nperm, 1) ||
      batch_size %% 1 != 0
  ) {
    stop("`batch_size` must be a whole number between 1 and `nperm`.")
  }

  batch_size <- min(batch_size, nperm)

  set.seed(seed) # let set.seed validate the seed
  set.seed(NULL)

  if (
    !is.vector(min_size, mode = "numeric") ||
      length(min_size) != 1L ||
      is.na(min_size) ||
      min_size < 2L ||
      min_size >= n_genes ||
      min_size %% 1 != 0
  ) {
    stop("`min_size` must be >= 2 and < nrow(X).")
  }

  if (
    !(is.infinite(max_size) || is.vector(max_size, mode = "numeric")) ||
      length(max_size) != 1L ||
      is.na(max_size) ||
      max_size < min_size ||
      (!is.infinite(max_size) && max_size %% 1 != 0)
  ) {
    stop("`max_size` must be Inf or an integer >= min_size.")
  }

  if (
    !is.vector(adjust_globally, mode = "logical") ||
      length(adjust_globally) != 1L ||
      is.na(adjust_globally)
  ) {
    stop("`adjust_globally` must be TRUE or FALSE.")
  }

  if (
    !is.vector(sort, mode = "logical") ||
      length(sort) != 1L ||
      is.na(sort)
  ) {
    stop("`sort` must be TRUE or FALSE.")
  }
}


#' @title Create Sparse Incidence Matrices
#'
#' @description Create a list of sparse incidence matrices, where the unique
#'   elements are rows and the unique sets are columns.
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
#' @importFrom data.table data.table := chmatch setattr
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
.sparseIncidence <- function(gene_sets, background) {
  if (
    !is.list(gene_sets) ||
      is.null(names(gene_sets))
  ) {
    stop("`gene_sets` must be a named list of character vectors.")
  }

  elements <- unlist(
    x = gene_sets,
    recursive = FALSE,
    use.names = FALSE
  )

  if (!is.vector(elements, mode = "character")) {
    stop("`gene_sets` must be a named list of character vectors.")
  }

  # "elements" is a factor, "sets" is not
  dt <- data.table(
    elements = elements,
    stringsAsFactors = TRUE
  )

  dt[, sets := rep.int(
    names(gene_sets),
    lengths(gene_sets)
  )]

  # Determine which elements are expected to be "down", if any
  dt[, direction_down := grepl(";d$", elements, perl = TRUE)]

  # Strip information about direction of change. This may reduce the number of
  # levels if an element is both "up" and "down": "gene;u" and "gene;d" become
  # "gene".
  setattr(dt$elements, "levels", sub(";[ud]{1}$", "", levels(dt[["elements"]])))

  # Do not chain with previous line, since the number of levels may change.
  unique_elements <- levels(dt[["elements"]])

  # Convert to characters to use chmatch()
  dt[, elements := as.character.factor(elements)]

  # Only need to check those elements of the background that overlap with the
  # elements of x. No need to validate background, since it gets validated
  # prior to calling this function.
  unique_elements <- intersect(unique_elements, background)

  if (length(unique_elements) == 0L) {
    stop("No elements of `gene_sets` are present in rownames(X).")
  }

  # Fast character matching (row indices for sparse matrix)
  dt[, i := chmatch(elements, unique_elements, nomatch = 0L)]

  if (any(dt[["i"]] == 0L)) {
    # Remove elements not in the background
    dt <- dt[i != 0L]

    unique_sets <- unique.default(dt[["sets"]])
  } else {
    unique_sets <- names(gene_sets)
  }

  # Column indices for sparse matrix
  dt[, j := chmatch(sets, unique_sets)]

  dim_names <- list(unique_elements, unique_sets)
  dims <- lengths(dim_names)

  # Keep genes expected to be "up"
  if (any(dt[["direction_down"]])) {
    dt_u <- dt[direction_down == FALSE]
  } else {
    dt_u <- dt
  }

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
  attr(A, which = "x") <- rep.int(1, length(attr(A, which = "x")))

  A_d <- NULL # default

  if (nrow(dt_u) < nrow(dt)) {
    dt_d <- dt[direction_down == TRUE]

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

    attr(A_d, which = "x") <- rep.int(1, length(attr(A_d, which = "x")))

    # The Hadamard product A * A.d should be a matrix of zeros
    if (length(attr(A * A_d, which = "x"))) {
      stop(
        "Elements can not be both up (suffix \";u\") and ",
        "down (suffix \";d\") in the same set."
      )
    }
  }

  out <- list(
    "A" = A,
    "A_d" = A_d
  )

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
.calcRankMatrix <- function(X) {
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
#' @param max_size integer or \code{Inf}; size of the largest gene set that will
#'   be tested.
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
#'   large gene sets removed.}
#'
#'   \item{"A_d"}{incidence matrix \code{A_d} with extremely small or extremely
#'   large gene sets removed.}
#' }
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.calcSetSize <- function(n,
                         Z_prime,
                         A,
                         A_d = NULL,
                         min_size = 2L,
                         max_size = Inf) {
  # Matrix with samples as rows and genes as columns. Elements are the number
  # of genes in each set with nonmissing values in the sample.
  M <- .Cpp_matmult_sparse(Z_prime, A) # Z'A
  M[M < min_size] <- 0L

  # It is extremely unlikely that a gene set would consist of all genes with
  # nonmissing values, but the set size limit will be 1 less than the smallest
  # number of nonmissing values.
  max_size <- max(min_size, min(max_size, min(n) - 1L))

  if (!is.null(A_d)) {
    M_d <- .Cpp_matmult_sparse(Z_prime, A_d)
    M_d[M_d < min_size] <- 0L

    small_sets_M <- apply(M == 0L, 2L, any)
    small_sets_M_d <- apply(M_d == 0L, 2L, any)

    # Combined set size exceeds max_size
    large_sets <- apply((M + M_d) > max_size, 2L, any)

    # Too small in both up and down portions or too large overall
    extreme_sets <- which((small_sets_M & small_sets_M_d) | large_sets)
  } else {
    M_d <- W_d <- NULL

    extreme_sets <- which(
      apply(M == 0 | M > max_size, 2L, any)
    )
  }

  if (length(extreme_sets)) {
    # If any sets are extreme, check if they are all extreme.
    if (length(extreme_sets) == ncol(A)) {
      stop(
        "All sets in `gene_sets` contain fewer than `min_size` genes or more ",
        "than `max_size` genes with nonmissing values."
      )
    }

    # Remove extreme sets
    A <- A[, -extreme_sets, drop = FALSE]
    M <- M[, -extreme_sets, drop = FALSE]

    if (!is.null(A_d)) {
      A_d <- A_d[, -extreme_sets, drop = FALSE]
      M_d <- M_d[, -extreme_sets, drop = FALSE]
    }
  }

  # Number of genes not in each set with nonmissing values
  W <- n - M

  if (!is.null(A_d)) {
    W_d <- n - M_d
  }

  out <- list(
    "A" = A,
    "A_d" = A_d,
    "M" = M,
    "W" = W,
    "M_d" = M_d,
    "W_d" = W_d
  )

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
#' @param A sparse incidence matrix with single genes as rows and gene sets as
#'   columns. Only those entries for genes expected to be up-regulated or
#'   lacking an expected direction of change are 1.
#' @param M matrix containing the number of nonmissing values in each set with
#'   samples as rows and gene sets as columns.
#' @param W matrix containing the number of nonmissing values not in each set
#'   with the same dimensions as \code{M}.
#' @param A_d,M_d,W_d like \code{A}, \code{M_d}, and \code{W_d}, but only for
#'   those genes expected to be down-regulated. If not testing directional sets,
#'   these will be \code{NULL}.
#'
#' @returns A list containing matrices of enrichment scores with names "ES",
#'   "ES_u", and "ES_d". Each matrix has \code{nrow(Y_prime)} samples as rows
#'   and \code{ncol(A)} gene sets as columns. If the sets are not directional,
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
                    W_d = NULL) {
  # Sample by gene set matrix of enrichment scores
  ES_u <- .Cpp_calcES(
    alpha,
    min_size,
    Y_prime,
    R_prime,
    sumRanks,
    A,
    M,
    W
  )

  if (!is.null(A_d)) { # directional database
    ES_d <- .Cpp_calcES(
      alpha,
      min_size,
      Y_prime,
      R_prime,
      sumRanks,
      A_d,
      M_d,
      W_d
    )

    ES <- ES_u - ES_d

    out <- list(
      "ES" = ES,
      "ES_u" = ES_u,
      "ES_d" = ES_d
    )
  } else {
    out <- list("ES" = ES_u)
  }

  return(out)
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
                            seed = NULL) {
  if (nperm != 0L) {
    n_batches <- ceiling(nperm / batch_size)

    batch_sizes <- c(
      rep.int(batch_size, n_batches - 1L),
      nperm - batch_size * (n_batches - 1L)
    )

    batch_id <- rep.int(seq_len(n_batches), batch_sizes)

    # Seeds for permutations
    set.seed(seed)
    seeds <- sample.int(1e7L, nperm, FALSE)

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
.szudzikPairing <- function(x, y) {
  # Do not need to validate x and y, since this function is only used with
  # integer vectors by other internal functions. Note that x * x prevents
  # integer overflow, unlike x ^ 2.
  x +
    ifelse(
      x < y,
      y * y, # y ^ 2 + x
      x * x + y # x ^ 2 + x + y
    )
}


#' @title Get Unique Set Sizes and Other Information for Permutations
#'
#' @description Construct a list containing vectors needed for permutations.
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
#'   This is used by \code{.makeResultsTable}.}
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
#'   \item{"theta_w_i"}{vector of the same length as \code{m_i}. Calculated as
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
.getUniqueSizes <- function(n_i,
                            m_i,
                            m_d_i = NULL) {
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

    # Each unique pair of entries in m_i and m_d_i are converted to a unique
    # integer. The same is done for each pair of entries in theta_m_i and
    # theta_m_d_i. Then, a vector is generated that maps each (m_i, m_d_i)
    # pair to the unique (theta_m_i, theta_m_d_i) pairs.
    rep_idx <- match(
      x = .szudzikPairing(m_i, m_d_i),
      table = .szudzikPairing(theta_m_i, theta_m_d_i)
    )
  } else {
    # Unique number of genes in each set
    theta_m_i <- unique(m_i)

    # Unique number of genes not in each set
    theta_w_i <- n_i - theta_m_i

    theta_m_d_i <- theta_w_d_i <- NULL

    max_set_size <- max(theta_m_i)

    rep_idx <- match(x = m_i, table = theta_m_i)
  }

  out <- list(
    "rep_idx" = rep_idx,
    "max_set_size" = max_set_size,
    "theta_m_i" = theta_m_i,
    "theta_w_i" = theta_w_i,
    # These may be NULL
    "theta_m_d_i" = theta_m_d_i,
    "theta_w_d_i" = theta_w_d_i
  )

  return(out)
}


.create_list <- function(n) {
  list(
    "n_same_sign" = rep.int(0L, n),
    "n_as_extreme" = rep.int(0L, n),
    "sum_ES_perm" = rep.int(0, n)
  )
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
#' @importFrom data.table data.table := setorderv rbindlist
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
                              ES_d_i) {
  # This will store the results for a single sample
  tab_i <- data.table(
    set = sets,
    set_size = m_i,
    m_d_i = m_d_i,
    ES = ES_i,
    row_order = seq_along(ES_i),
    stringsAsFactors = FALSE
  )

  # Sorting by set size and then ES allows for major optimizations
  setorderv(
    x = tab_i,
    cols = c("set_size", "ES", "row_order"),
    order = c(1L, 1L, 1L)
  )

  # Reordered
  m_i <- tab_i[["set_size"]]
  m_d_i <- tab_i[["m_d_i"]]

  if (nperm == 0L) {
    tab_i[, `:=`(
      n_same_sign = rep.int(0L, length(ES_i)),
      n_as_extreme = rep.int(0L, length(ES_i)),
      sum_ES_perm = rep.int(0, length(ES_i))
    )]
  } else {
    # Unique set sizes and other information for permutations
    size_list <- .getUniqueSizes(
      n_i = n_i,
      m_i = m_i,
      m_d_i = m_d_i
    )

    # Extract list components: rep_idx, max_set_size, theta_m_i, theta_w_i,
    # theta_m_d_i, and theta_w_d_i
    list2env(x = size_list, envir = environment())

    ES_ls <- split(
      x = tab_i[["ES"]],
      f = rep_idx
    )

    names(ES_ls) <- NULL

    # Indices of non-missing values for a particular column
    element_indices <- which(r_i != 0L)

    if (length(element_indices) != length(r_i)) {
      r_i <- r_i[element_indices]
      y_i <- y_i[element_indices]
    }

    # List to store results needed to calculate P-values and NES
    perm_ls <- lapply(table(rep_idx), .create_list)
    names(perm_ls) <- NULL

    if (is.null(theta_m_d_i)) {
      # Split permutations into batches to reduce memory consumption
      for (seeds in seed_list) {
        ES_perm <- .Cpp_calcESPerm(
          alpha,
          y_i,
          r_i,
          seeds,
          max_set_size,
          sumRanks_i,
          theta_m_i,
          theta_w_i
        )

        # Update perm_ls by reference
        .Cpp_extractPermInfo(perm_ls, ES_ls, ES_perm)
      }
    } else {
      for (seeds in seed_list) {
        # Directional gene sets
        ES_perm <- .Cpp_calcESPerm_dir(
          alpha,
          y_i,
          r_i,
          seeds,
          max_set_size,
          sumRanks_i,
          theta_m_i,
          theta_w_i,
          theta_m_d_i,
          theta_w_d_i,
          min_size
        )

        # Update perm_ls by reference
        .Cpp_extractPermInfo(perm_ls, ES_ls, ES_perm)
      }
    } # end permutation batching

    # Convert each list of vectors to a data.table to use rbindlist
    for (i in seq_len(length(perm_ls))) {
      class(perm_ls[[i]]) <- "data.table"
    }

    perm_dt <- rbindlist(perm_ls, use.names = FALSE)

    # Make columns for summary vectors in results
    tab_i[, `:=`(
      n_same_sign = perm_dt[["n_same_sign"]],
      n_as_extreme = perm_dt[["n_as_extreme"]],
      sum_ES_perm = perm_dt[["sum_ES_perm"]]
    )]
  } # end permutations

  setorderv(tab_i, cols = "row_order", order = 1L) # original row order

  if (!is.null(m_d_i)) {
    tab_i[, `:=`(
      set_size = set_size + m_d_i,
      ES_u = ES_u_i,
      ES_d = ES_d_i
    )]
  }

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
#' @param alternative character; the alternative hypothesis. One of
#'   "\code{two.sided}" (default), "\code{less}", or "\code{greater}". The
#'   latter two will perform one-sided tests.
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
#' @importFrom data.table rbindlist := setorderv setDF
#' @importFrom stats p.adjust
#'
#' @noRd
.stackResults <- function(tab,
                          nperm = 1000L,
                          sort = TRUE,
                          adjust_globally = FALSE,
                          alternative = "two.sided") {
  sample_names <- names(tab)
  tab <- rbindlist(tab, use.names = FALSE, idcol = "sample")

  tab[, `:=`(
    sample = factor(sample, levels = sample_names),
    set_size = as.integer(set_size),
    NES = ES / (sum_ES_perm / n_same_sign)
  )]

  switch(
    EXPR = alternative,
    two.sided = {
      tab[, p_value := (n_as_extreme + 1L) / (n_same_sign + 1L)]
    },
    less = {
      tab[, p_value := ifelse(
        ES < 0,
        n_as_extreme + 1L,
        nperm - n_as_extreme + 1L
      ) / (nperm + 1L)]
    },
    greater = {
      tab[, p_value := ifelse(
        ES >= 0,
        n_as_extreme + 1L,
        nperm - n_as_extreme + 1L
      ) / (nperm + 1L)]
    }
  )

  tab[
    n_same_sign == 0L, # only happens when nperm is very small
    `:=`(
      NES = NA_real_,
      p_value = NA_real_
    )
  ]

  if (nperm == 0L) {
    tab[, `:=`(
      n_same_sign = NA_integer_,
      n_as_extreme = NA_integer_,
      p_value = NA_real_
    )]
  } else if (sort) {
    setorderv(x = tab, cols = c("sample", "p_value"), order = c(1L, 1L))
  }

  if (adjust_globally) {
    tab[, adj_p_value := p.adjust(p_value, method = "BH")]
  } else {
    tab[, adj_p_value := p.adjust(p_value, method = "BH"), by = sample]
  }

  # Reorder/select columns
  keep_cols <- intersect(
    c(
      "sample", "set", "set_size",
      "ES_u", "ES_d", "ES", "NES",
      "n_same_sign", "n_as_extreme",
      "p_value", "adj_p_value"
    ),
    colnames(tab)
  )

  tab <- tab[, keep_cols, with = FALSE]

  # Convert to data.frame
  setDF(tab)

  return(tab)
}
