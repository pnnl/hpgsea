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
                            nperm = 1e5L,
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
      nperm > 2e9L || # arbitrary limit on number of permutations
      nperm %% 1 != 0 # decimal number
  ) {
    stop("`nperm` must be a whole number between 0 and 2 billion.")
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
#' @param min_size integer; minimum gene set size required for testing. Default
#'   is 2.
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
#' @importFrom collapse alloc allv anyv fsubset funique groupid vec vlengths
#'   vtypes whichNA whichv
#' @importFrom data.table chmatch
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
.sparseIncidence <- function(gene_sets,
                             background,
                             min_size = 2L) {
  # background and min_size are validated before calling this function

  err <- "`gene_sets` must be a named list of character vectors."

  if (!is.list(gene_sets)) {
    stop(err)
  }

  if (is.null(names(gene_sets))) {
    stop(err)
  }

  all_char <- allv(vtypes(gene_sets, use.names = FALSE), "character")

  if (!all_char) {
    stop(err)
  }

  # Pre-filter to remove gene sets that are too small. We can not remove gene
  # sets that are too large without first restricting the genes to the
  # background. Additionally, missing values in matrix X may make the gene sets
  # smaller than min_size, and duplicate elements may cause gene sets to contain
  # at least min_size genes. In all of these cases, the gene sets will survive
  # the initial filter, but they will be removed later by .calcSetSize.
  set_sizes <- vlengths(gene_sets)

  keep_sets <- whichv(set_sizes >= min_size, TRUE)

  if (length(keep_sets) != length(gene_sets)) {
    if (length(keep_sets) == 0L) {
      stop(
        "No gene sets with at least `min_size` elements."
      )
    }

    gene_sets <- gene_sets[keep_sets]
    set_sizes <- fsubset(set_sizes, keep_sets)
  }

  elements <- vec(gene_sets)
  unique_elements <- funique(elements)

  # Determine if any elements have an expected direction of change (i.e., they
  # end in ";u" or ";d")
  any_dir <- anyv(grepl(";[ud]{1}", unique_elements, perl = TRUE), TRUE)

  if (any_dir) {
    # Convert elements to an integer vector to index direction_down
    elements <- chmatch(elements, unique_elements)

    # Determine which elements are expected to be "down", if any. grepl will be
    # slower than indexing, so apply grepl to the unique elements and then
    # expand the logical vector by indexing with the elements integer vector.
    direction_down <- grepl(";d$", unique_elements, perl = TRUE)
    direction_down <- fsubset(direction_down, elements)

    # Strip information about direction of change. This may reduce the number of
    # levels if an element is both "up" and "down": "gene;u" and "gene;d" become
    # "gene".
    unique_elements <- sub(";[ud]{1}$", "", unique_elements)

    # Convert back to a character vector to use chmatch()
    elements <- fsubset(unique_elements, elements)
  } else {
    direction_down <- NULL # signal that this should not be used later
  }

  # Only need to check those elements of the background that overlap with the
  # elements of x.
  unique_elements <- intersect(background, unique_elements)

  if (length(unique_elements) == 0L) {
    stop("No elements of `gene_sets` are present in rownames(X).")
  }

  # Row indices for sparse matrix
  i <- chmatch(elements, unique_elements, nomatch = NA_integer_)

  # Unique set names
  unique_sets <- names(gene_sets)

  # Column indices for sparse matrix
  j <- rep.int(seq_along(unique_sets), set_sizes)

  if (anyNA(i)) {
    # Remove elements not in the background
    keep <- whichNA(i, invert = TRUE)

    i <- fsubset(i, keep)
    j <- fsubset(j, keep)

    unique_sets <- fsubset(unique_sets, funique(j))

    j <- groupid(j) # this works because each set is a contiguous batch
  }

  dim_names <- list(unique_elements, unique_sets)
  dims <- vlengths(dim_names)

  # Split elements by direction of change to create two incidence matrices
  if (any_dir) {
    idx_down <- whichv(direction_down, TRUE)
    direction_down <- NULL # signal that this is no longer needed

    i_down <- fsubset(i, idx_down)
    j_down <- fsubset(j, idx_down)

    i <- fsubset(i, -idx_down)
    j <- fsubset(j, -idx_down)
  }

  # Incidence matrix where a 1 indicates that the element is in the set. If
  # gene_sets is a directional database, then A will only contain elements that
  # are expected to be "up".
  A <- sparseMatrix(
    i = i,
    j = j,
    x = alloc(1, length(i)),
    dims = dims,
    dimnames = dim_names,
    check = FALSE,
    use.last.ij = FALSE
  )

  # In the unlikely event where an element appears multiple times in the same
  # set, some values of A will be > 1. Replace all values with 1. Could also use
  # the use.last.ij parameter in sparseMatrix(), but this is faster.
  attr(A, which = "x") <- alloc(1, length(attr(A, which = "x")))

  A_d <- NULL # default

  if (any_dir) {
    # Incidence matrix where a 1 indicates that an element is expected to be
    # down in the set
    A_d <- sparseMatrix(
      i = i_down,
      j = j_down,
      x = alloc(1, length(i_down)),
      dims = dims,
      dimnames = dim_names,
      check = FALSE,
      use.last.ij = FALSE
    )

    attr(A_d, which = "x") <- alloc(1, length(attr(A_d, which = "x")))

    # The Hadamard product A * A_d should be a matrix of zeros, since genes can
    # not be "up" and "down" in the same set.
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
.calcES <- function(min_size = 2L,
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
#' @param n the total number of nonmissing values in the i-th sample.
#' @param m integer vector containing the number of genes in each set or the
#'   number of genes expected to be up-regulated in each set in the i-th sample.
#' @param m_d integer vector or \code{NULL}. The number of genes expected to be
#'   down-regulated in each set in the i-th sample. This vector has the same
#'   length as \code{m}.
#'
#' @returns A named list.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.getUniqueSizes <- function(n,
                            m,
                            m_d = NULL) {
  if (is.null(m_d)) {
    # Level 1 ----

    # All gene sets, even if the number of genes is not unique (m). Set sizes
    # are sorted.

    # Level 2 ----

    # Unique number of genes in each set
    L2_m <- unique(m)

    # Unique number of genes not in each set
    L2_w <- n - L2_m

    # Indices to map from level 2 to level 1 ----
    map_L2_to_L1 <- match(m, L2_m)

    out <- list(
      "max_set_size" = max(L2_m),
      "L2_m" = L2_m,
      "L2_w" = L2_w
    )
  } else {
    # Level 1 ----

    # All gene sets, even if the number of up and down-regulated genes in the
    # set is not unique (m and m_d).

    # Level 2 ----

    # Unique combinations of the number of up and down genes. Sorted by the
    # number of up-regulated genes and then by the number of down-regulated
    # genes.
    unique_pairs <- unique(cbind(m, m_d))

    # Up-regulated, in the set
    L2_m <- unique_pairs[, 1L]

    # Not up-regulated (includes genes not in the set)
    L2_w <- theta_w <- n - L2_m

    # Down-regulated, in the set
    L2_m_d <- unique_pairs[, 2L]

    # Not down-regulated (includes genes not in the set)
    L2_w_d <- n - L2_m_d

    # Level 3 ----

    # Unique number of up-regulated genes and the unique number of
    # down-regulated genes (separate, sorted vectors). Same as level 2 for
    # non-directional gene sets.

    # Up-regulated, in the set
    L3_m <- unique(L2_m)

    # Not up-regulated (includes genes not in the set)
    L3_w <- n - L3_m

    # Down-regulated, in the set
    L3_m_d <- sort(unique(L2_m_d))

    # Not down-regulated (includes genes not in the set)
    L3_w_d <- n - L3_m_d

    # Indices to map from lower levels to higher levels ----
    map_L2_to_L1 <- match(
      .szudzikPairing(m, m_d),
      .szudzikPairing(L2_m, L2_m_d)
    )

    map_L3_to_L2 <- match(L2_m, L3_m) # up-regulated
    map_L3_to_L2_d <- match(L2_m_d, L3_m_d) # down-regulated

    out <- list(
      "max_set_size" = max(L2_m + L2_m_d),
      "L3_m" = L3_m,
      "L3_w" = L3_w,
      "L3_m_d" = L3_m_d,
      "L3_w_d" = L3_w_d,
      "map_L3_to_L2" = map_L3_to_L2,
      "map_L3_to_L2_d" = map_L3_to_L2_d
    )
  }

  out[["ES_end"]] <- which(!duplicated(map_L2_to_L1, fromLast = TRUE))

  return(out)
}


# TODO update documentation

#' @title Generate ssGSEA Results Table for a Single Sample
#'
#' @inheritParams fast_ssgsea
#' @param y numeric vector of absolute gene-level values from the i-th sample
#'   raised to the power of \code{alpha}.
#' @param r numeric vector of ranks of the gene-level values from the i-th
#'   sample.
#' @param n integer; number of genes in the i-th sample with nonmissing values.
#' @param sumRanks integer; the sum of \code{r}.
#' @param m the number of genes with nonmissing values in each set from the i-th
#'   sample.
#' @param m_d the number of genes with nonmissing values not in each set from
#'   the i-th sample.
#' @param sets character vector of gene set labels.
#' @param ES numeric vector of enrichment scores from the i-th sample.
#' @param ES_u numeric vector of enrichment scores for the genes in each set
#'   that are expected to be up-regulated in the i-th sample. If not testing a
#'   directional database, this will be \code{NULL}.
#' @param ES_d numeric vector of enrichment scores for the genes in each set
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
.makeResultsTable <- function(seed = NULL,
                              nperm = 1e5L,
                              batch_size = 1000L,
                              min_size = 2L,
                              y,
                              r,
                              n,
                              sumRanks,
                              m,
                              m_d,
                              sets,
                              ES,
                              ES_u,
                              ES_d) {
  tab_i <- data.table(
    set = sets,
    set_size = m,
    m_d = m_d,
    ES = ES,
    row_order = seq_along(ES),
    stringsAsFactors = FALSE
  )

  # Sorting by set size and then ES allows for major optimizations
  setorderv(
    x = tab_i,
    cols = intersect(
      c("set_size", "m_d", "ES", "row_order"),
      colnames(tab_i)
    )
  )

  # Reordered
  m <- tab_i[["set_size"]]
  m_d <- tab_i[["m_d"]]
  ES <- tab_i[["ES"]]

  # Vectors needed to store information from the permutation tests to calculate
  # NES and p-values
  n_same_sign <- rep.int(0L, length(ES))
  n_as_extreme <- rep.int(0L, length(ES))
  sum_ES_perm <- rep.int(0, length(ES))

  if (nperm != 0L) {
    size_list <- .getUniqueSizes(
      n = n,
      m = m,
      m_d = m_d
    )

    list2env(x = size_list, envir = environment())

    # Indices of non-missing values for a particular column
    element_indices <- which(r != 0)

    if (length(element_indices) != length(r)) {
      r <- r[element_indices]
      y <- y[element_indices]
    }

    if (is.null(m_d)) {
      .Cpp_calc_ES_perm(
        n_same_sign,
        n_as_extreme,
        sum_ES_perm,
        seed,
        nperm,
        batch_size,
        ES,
        ES_end,
        y,
        r,
        max_set_size,
        sumRanks,
        L2_m,
        L2_w
      )
    } else {
      .Cpp_calc_ES_perm_dir(
        n_same_sign,
        n_as_extreme,
        sum_ES_perm,
        seed,
        nperm,
        batch_size,
        ES,
        ES_end,
        y,
        r,
        max_set_size,
        sumRanks,
        L3_m,
        L3_w,
        L3_m_d,
        L3_w_d,
        map_L3_to_L2,
        map_L3_to_L2_d
      )
    }
  } # end permutations

  tab_i[, `:=`(
    n_same_sign = n_same_sign,
    n_as_extreme = n_as_extreme,
    sum_ES_perm = sum_ES_perm
  )]

  setorderv(tab_i, cols = "row_order", order = 1L) # original row order

  if (!is.null(m_d)) {
    tab_i[, `:=`(
      set_size = set_size + m_d,
      ES_u = ES_u,
      ES_d = ES_d
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
                          nperm = 1e5L,
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

  # NOTE this shouldn't happen. Throw an error instead?
  tab[
    is.na(NES) | is.infinite(NES),
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
