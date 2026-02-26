#' @title Validate Gene-Level Values
#'
#' @description Validate `stats`, remove missing values, and sort genes
#'   lexicographically.
#'
#' @inheritParams fast_ssgsea
#'
#' @author Tyler Sagendorf
#'
#' @importFrom collapse any_duplicated radixorderv whichNA
#'
#' @noRd
.prepare_stats <- function(stats) {
  if (
    !is.vector(stats, mode = "numeric") ||
      is.null(names(stats)) ||
      anyNA(names(stats)) ||
      any_duplicated(names(stats))
  ) {
    stop("`stats` must be a numeric vector with unique names.", call. = FALSE)
  }

  # Remove missing values
  idx_not_NA <- whichNA(stats, invert = TRUE)

  if (length(idx_not_NA) < 3L) {
    stop("`stats` must have at least 3 nonmissing values.", call. = FALSE)
  }

  genes <- fsubset(names(stats), idx_not_NA)
  stats <- fsubset(stats, idx_not_NA)

  # Sort genes lexicographically so results are consistent when there are ties
  o <- radixorderv(genes)
  stats <- fsubset(stats, o)
  names(stats) <- fsubset(genes, o)

  storage.mode(stats) <- "double"

  return(stats)
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
.validate_params <- function(alpha = 1,
                             nperm = 1e5L,
                             min_size = 2L,
                             max_size = Inf,
                             sort = TRUE,
                             seed = NULL,
                             n_genes) {
  if (
    !is.vector(alpha, mode = "numeric") ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      is.infinite(alpha) ||
      alpha < 0
  ) {
    stop("`alpha` must be a single non-negative real number.", call. = FALSE)
  }

  if (
    !is.vector(nperm, mode = "numeric") ||
      length(nperm) != 1L ||
      is.na(nperm) ||
      nperm < 0L ||
      nperm > 2e9L || # close to .Machine$integer.max
      nperm %% 1 != 0
  ) {
    stop("`nperm` must be an integer between 0 and 2 billion.", call. = FALSE)
  }

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
    stop(
      "`min_size` must be >= 2 and less than ",
      "the number of nonmissing values in `stats`.",
      call. = FALSE
    )
  }

  if (
    !(is.infinite(max_size) || is.vector(max_size, mode = "numeric")) ||
      length(max_size) != 1L ||
      is.na(max_size) ||
      max_size < min_size ||
      (!is.infinite(max_size) && max_size %% 1 != 0)
  ) {
    stop("`max_size` must be Inf or an integer >= min_size.", call. = FALSE)
  }

  if (
    !is.vector(sort, mode = "logical") ||
      length(sort) != 1L ||
      is.na(sort)
  ) {
    stop("`sort` must be TRUE or FALSE.", call. = FALSE)
  }
}


#' @title Fast, specialized rep.int
#'
#' @description Equivalent to `rep.int(seq_along(times), times)`, but several
#' times faster for when the output is large.
#'
#' @param times integer vector of group sizes.
#'
#' @returns Integer vector.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.C_rep_int <- function(sizes) {
  .Call("_C_rep_int", sizes)
}


#' @title Calculate group sizes
#'
#' @description Given a vector of integers from 1 to at most `n_groups`, count
#'   the number of occurrences in each group. Like `table`, but the counts can
#'   be 0.
#'
#' @param group_ids vector of integers from 1 to at most `n_groups`.
#' @param n_groups integer; number of groups. At least `max(group_ids)`, but
#'   this is not checked.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.C_group_sizes <- function(group_ids, n_groups) {
  .Call("_C_group_sizes", group_ids, n_groups)
}


#' @title Calculate Enrichment Scores
#'
#' @param y_prime numeric vector of absolute gene-level values raised to the
#'   power of `alpha` for genes that are members of at least one gene set.
#' @param r_prime numeric vector of the ranks of the gene-level values for genes
#'   that are members of at least one gene set.
#' @param sum_ranks numeric; the sum of all ranks.
#' @param i integer vector; indices of the genes in all sets. Used to index
#'   vectors `y_prime` and `r_prime`.
#' @param m integer vector; the number of genes in each set. Used to select
#'   elements of `i`.
#' @param w integer vector; the number of genes that are not in each set.
#' @inheritParams fast_ssgsea
#'
#' @returns Numeric vector of enrichment scores with the same length as `m`. If
#'   all elements of `y_prime` are 0 for a particular set, the ES for that set
#'   will be `NA`. If any `m` are less than `min_size`, the correspond ES will
#'   be 0.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.C_calc_ES <- function(y_prime, r_prime, sum_ranks, i, m, w, min_size) {
  .Call("_C_calc_ES", y_prime, r_prime, sum_ranks, i, m, w, min_size)
}


#' @title Calculate Enrichment Scores and Other Information
#'
#' @description Calculate enrichment scores and information needed for
#'   permutation tests later on.
#'
#' @inheritParams fast_ssgsea
#' @param n_genes integer; total number of genes. The length of `stats`.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom collapse %!iin% allv anyv fmatch fmax fsubset funique groupid vec
#'   vlengths vtypes whichNA whichv
#' @importFrom data.table frank
#'
#' @noRd
.calc_ES <- function(stats,
                     alpha = 1,
                     n_genes,
                     gene_sets,
                     min_size = 2L,
                     max_size = Inf) {
  max_size <- max(min_size, min(max_size, n_genes - 1L))

  storage.mode(min_size) <- storage.mode(max_size) <- "integer"

  err <- "`gene_sets` must be a named list of character vectors."

  if (!is.list(gene_sets)) {
    stop(err, call. = FALSE)
  }

  if (is.null(names(gene_sets))) {
    stop(err, call. = FALSE)
  }

  all_char <- allv(vtypes(gene_sets, use.names = FALSE), "character")

  if (!all_char) {
    stop(err, call. = FALSE)
  }

  # Pre-filter to remove gene sets that are too small. We can not remove gene
  # sets that are too large without first restricting the genes to the
  # background. Additionally, duplicate elements may cause gene sets to contain
  # at least min_size genes. In this case, the gene sets will survive the
  # initial filter, but they will be removed toward the end of this function.
  set_sizes <- vlengths(gene_sets)

  keep_sets <- whichv(set_sizes >= min_size, TRUE)

  if (length(keep_sets) != length(gene_sets)) {
    if (length(keep_sets) == 0L) {
      stop("No gene sets with at least `min_size` elements.", call. = FALSE)
    }

    gene_sets <- gene_sets[keep_sets]
    set_sizes <- fsubset(set_sizes, keep_sets)
  }

  elements <- vec(gene_sets)
  unique_elements <- funique(elements)

  # Determine if any elements have an expected direction of change
  any_dir <- anyv(grepl(";[ud]{1}", unique_elements, perl = TRUE), TRUE)

  if (any_dir) {
    # Convert elements to an integer vector to index direction_down
    elements <- fmatch(elements, unique_elements)

    # Determine which elements are expected to be "down", if any. grepl will be
    # slower than indexing, so apply grepl to the unique elements and then
    # expand the logical vector by indexing with the elements integer vector.
    direction_down <- grepl(";d$", unique_elements, perl = TRUE)
    direction_down <- fsubset(direction_down, elements)

    # Strip information about direction of change. This may reduce the number of
    # levels if an element is both "up" and "down": "gene;u" and "gene;d" become
    # "gene".
    unique_elements <- sub(";[ud]{1}$", "", unique_elements)

    # Convert elements back to a character vector for fmatch
    elements <- fsubset(unique_elements, elements)
  } else {
    direction_down <- NULL # signal that this should not be used later
  }

  # Only need to check those elements of the background that overlap with the
  # elements of x
  unique_elements <- intersect(names(stats), unique_elements)

  if (length(unique_elements) == 0L) {
    stop(
      "No elements of `gene_sets` are present in names(stats).",
      call. = FALSE
    )
  }

  # Indices of the genes in each set
  i <- fmatch(elements, unique_elements, nomatch = NA_integer_)

  unique_sets <- names(gene_sets)

  # Unique integers for each set
  j <- .C_rep_int(set_sizes)

  if (anyNA(i)) {
    # Remove elements not in the background
    keep <- whichNA(i, invert = TRUE)

    i <- fsubset(i, keep)
    j <- fsubset(j, keep)

    if (any_dir) {
      direction_down <- fsubset(direction_down, keep)
    }

    unique_sets <- fsubset(unique_sets, funique(j))

    j <- groupid(j)
    class(j) <- "integer" # necessary for .C_group_sizes()
  }

  # Split by direction of change
  if (any_dir) {
    idx_down <- whichv(direction_down, TRUE)
    direction_down <- NULL # signal that this is no longer needed

    if (length(idx_down)) {
      i_down <- fsubset(i, idx_down)
      j_down <- fsubset(j, idx_down)

      i <- fsubset(i, -idx_down)
      j <- fsubset(j, -idx_down)
    } else {
      i_down <- j_down <- integer(0L)
    }
  }

  n_sets <- length(unique_sets)

  # Number of genes in each set, or the number of genes expected to be
  # up-regulated in each set
  m <- .C_group_sizes(j, n_sets)

  if (any_dir) {
    # Number of genes expected to be down-regulated in each set
    m_d <- .C_group_sizes(j_down, n_sets)

    extreme_sets <- which(
      (m < min_size & m_d < min_size) | (m + m_d) > max_size
    )
  } else {
    m_d <- NULL

    extreme_sets <- which(
      m < min_size | m > max_size
    )
  }

  # Remove gene sets that are too small or too large
  if (length(extreme_sets)) {
    if (length(extreme_sets) == n_sets) {
      stop(
        "All sets in `gene_sets` have fewer than `min_size` ",
        "or more than `max_size` genes in `stats`.",
        call. = FALSE
      )
    }

    unique_sets <- fsubset(unique_sets, -extreme_sets)
    m <- fsubset(m, -extreme_sets)
    i <- fsubset(i, j %!iin% extreme_sets)

    if (any_dir) {
      m_d <- fsubset(m_d, -extreme_sets)
      i_down <- fsubset(i_down, j_down %!iin% extreme_sets)
    }
  }

  # The number of genes not in each set, or the number of genes that are not in
  # the set and expected to be up-regulated
  w <- n_genes - m

  if (any_dir) {
    # A smaller max_size will speed up the permutation tests
    max_size <- fmax(m + m_d)

    # The number of genes that are not in the set and expected to be
    # down-regulated
    w_d <- n_genes - m_d
  } else {
    max_size <- fmax(m)

    w_d <- NULL
  }

  # The sum of all ranks is the n-th triangular number
  sum_ranks <- n_genes / 2 * (n_genes + 1)

  y <- abs(stats)^alpha
  r <- frank(stats, ties.method = "last")
  names(r) <- names(y)
  storage.mode(r) <- "double"

  # Subset to the genes that are elements of at least one gene set
  y_prime <- y[unique_elements]
  r_prime <- r[unique_elements]

  if (any_dir) {
    # Calculate enrichment scores separately for the up-regulated and
    # down-regulated genes
    ES_u <- .C_calc_ES(
      y_prime = y_prime,
      r_prime = r_prime,
      sum_ranks = sum_ranks,
      i = i,
      m = m,
      w = w,
      min_size = min_size
    )

    ES_d <- .C_calc_ES(
      y_prime = y_prime,
      r_prime = r_prime,
      sum_ranks = sum_ranks,
      i = i_down,
      m = m_d,
      w = w_d,
      min_size = min_size
    )

    ES <- ES_u - ES_d
  } else {
    ES <- .C_calc_ES(
      y_prime = y_prime,
      r_prime = r_prime,
      sum_ranks = sum_ranks,
      i = i,
      m = m,
      w = w,
      min_size = min_size
    )

    ES_u <- ES_d <- NULL
  }

  out <- list(
    "y" = y,
    "r" = r,
    "min_size" = min_size,
    "max_size" = max_size,
    "sets" = unique_sets,
    "sum_ranks" = sum_ranks,
    "ES" = ES,
    "ES_u" = ES_u,
    "ES_d" = ES_d,
    "m" = m,
    "w" = w,
    "m_d" = m_d,
    "w_d" = w_d
  )

  return(out)
}


#' @title Szudzik's Pairing Function for Two Integer Vectors
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
.pair_szudzik <- function(x, y) {
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
#' @inheritParams .calc_ES_perm
#'
#' @returns A named list.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom collapse fmatch funique
#'
#' @noRd
.unique_set_sizes <- function(n_genes,
                              m,
                              m_d = NULL) {
  if (is.null(m_d)) {
    # Level 1 ----

    # All gene sets, even if the number of genes is not unique (m). Set sizes
    # are sorted.

    # Level 2 ----

    # Unique number of genes in each set
    L2_m <- funique(m)

    # Unique number of genes not in each set
    L2_w <- n_genes - L2_m

    # Indices to map from level 2 to level 1 ----
    map_L2_to_L1 <- fmatch(m, L2_m)

    out <- list(
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
    L2_w <- n_genes - L2_m

    # Down-regulated, in the set
    L2_m_d <- unique_pairs[, 2L]

    # Not down-regulated (includes genes not in the set)
    L2_w_d <- n_genes - L2_m_d

    # Level 3 ----

    # Unique number of up-regulated genes and the unique number of
    # down-regulated genes (separate, sorted vectors). Same as level 2 for
    # non-directional gene sets.

    # Up-regulated, in the set
    L3_m <- funique(L2_m)

    # Not up-regulated (includes genes not in the set)
    L3_w <- n_genes - L3_m

    # Down-regulated, in the set
    L3_m_d <- funique(L2_m_d, sort = TRUE)

    # Not down-regulated (includes genes not in the set)
    L3_w_d <- n_genes - L3_m_d

    # Indices to map from lower levels to higher levels ----
    map_L2_to_L1 <- fmatch(
      .pair_szudzik(m, m_d),
      .pair_szudzik(L2_m, L2_m_d)
    )

    map_L3_to_L2 <- fmatch(L2_m, L3_m) # up-regulated
    map_L3_to_L2_d <- fmatch(L2_m_d, L3_m_d) # down-regulated

    out <- list(
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


#' @title Calculate Permutation Enrichment Scores
#'
#' @inheritParams fast_ssgsea
#' @inheritParams .calc_ES
#' @param ES_list list; output of `.calc_ES()`.
#'
#' @returns A `data.table`.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom collapse alloc fsubset whichNA
#' @importFrom data.table data.table := setorderv
#'
#' @noRd
.calc_ES_perm <- function(seed = NULL,
                          nperm = 1e5L,
                          n_genes,
                          ES_list) {
  list2env(ES_list, envir = environment())

  tab <- data.table(
    set = sets,
    m = m,
    m_d = m_d,
    ES = ES,
    ES_u = ES_u,
    ES_d = ES_d,
    row_order = seq_along(ES),
    stringsAsFactors = FALSE
  )

  if (is.null(m_d)) {
    tab[, set_size := m]
  } else {
    tab[, set_size := m + m_d]
  }

  # Sorting by set size and then ES allows for major optimizations
  setorderv(
    x = tab,
    cols = intersect(
      c("m", "m_d", "ES", "row_order"),
      colnames(tab)
    )
  )

  # Reordered
  m <- tab[["m"]]
  m_d <- tab[["m_d"]]
  ES <- tab[["ES"]]

  # The ES will be NA if the y for all genes in the set are 0
  idx_not_NA <- whichNA(ES, invert = TRUE)

  if (length(idx_not_NA) != length(ES)) {
    if (length(idx_not_NA) == 0L) {
      nperm <- 0L
    } else {
      m <- fsubset(m, idx_not_NA)
      ES <- fsubset(ES, idx_not_NA)

      if (!is.null(m_d)) {
        m_d <- fsubset(m_d, idx_not_NA)
      }
    }
  }

  # Vectors needed to store information from the permutation tests to calculate
  # NES and p-values
  n_same_sign <- alloc(0L, length(ES))
  n_as_extreme <- alloc(0L, length(ES))
  sum_ES_perm <- alloc(0, length(ES))

  if (nperm != 0L) {
    size_list <- .unique_set_sizes(
      n_genes = n_genes,
      m = m,
      m_d = m_d
    )

    list2env(x = size_list, envir = environment())

    if (is.null(m_d)) {
      .Cpp_calc_ES_perm(
        n_same_sign,
        n_as_extreme,
        sum_ES_perm,
        seed,
        nperm,
        ES,
        ES_end,
        y,
        r,
        max_size,
        sum_ranks,
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
        ES,
        ES_end,
        y,
        r,
        max_size,
        sum_ranks,
        L3_m,
        L3_w,
        L3_m_d,
        L3_w_d,
        map_L3_to_L2,
        map_L3_to_L2_d
      )
    }
  }

  tab[idx_not_NA, `:=`(
    n_same_sign = n_same_sign,
    n_as_extreme = n_as_extreme,
    sum_ES_perm = sum_ES_perm
  )]

  setorderv(tab, cols = "row_order", order = 1L) # original row order

  return(tab)
}


#' @title Calculate P-values and Normalized Enrichment Scores
#'
#' @param tab a `data.table`. The output of `.calc_ES_perm()`.
#' @param nperm integer; the number of permutations.
#' @param sort logical; whether to sort rows in descending order by p-value.
#' @param alternative character; the alternative hypothesis. One of "two.sided"
#'   (default), "less", or "greater". The latter two will perform one-sided
#'   tests.
#'
#' @returns A `data.frame`.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom collapse whichNA
#' @importFrom data.table := setorderv setDF
#' @importFrom stats p.adjust
#'
#' @noRd
.calc_pvals <- function(tab,
                        nperm = 1e5L,
                        sort = TRUE,
                        alternative = "two.sided") {
  tab[, `:=`(
    set_size = set_size,
    NES = ES / (sum_ES_perm / n_same_sign)
  )]

  switch(
    EXPR = alternative,
    two.sided = {
      tab[
        whichNA(ES, invert = TRUE),
        p_value := (n_as_extreme + 1L) / (n_same_sign + 1L)
      ]
    },
    less = {
      tab[whichNA(ES, invert = TRUE), p_value := ifelse(
        ES < 0,
        n_as_extreme + 1L,
        nperm - n_as_extreme + 1L
      ) / (nperm + 1L)]
    },
    greater = {
      tab[whichNA(ES, invert = TRUE), p_value := ifelse(
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
    setorderv(x = tab, cols = "p_value", order = 1L)
  }

  tab[, adj_p_value := p.adjust(p_value, method = "BH")]

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

  setDF(tab)

  return(tab)
}
