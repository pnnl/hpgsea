#' @title Validate Gene-Level Values
#'
#' @description Validate `stats`, remove missing values, and sort genes
#'   lexicographically.
#'
#' @inheritParams hpgsea
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


#' @title Validate hpgsea function parameters
#'
#' @inheritParams hpgsea
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


#' @title Convert gene sets to a list of indices for `.calc_ES`
#'
#' @inheritParams hpgsea
#'
#' @returns A named list.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom collapse allv anyv fmatch fsubset funique groupid vec vlengths
#' @importFrom collapse vtypes whichNA whichv
#'
#' @noRd
.gene_sets_to_indices <- function(stats,
                                  gene_sets,
                                  min_size = 2L) {
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
  # sets that are too large without first restricting genes to names(stats).
  set_sizes <- vlengths(gene_sets)

  keep_sets <- whichv(set_sizes >= min_size, TRUE)

  if (length(keep_sets) == 0L) {
    stop("No gene sets with at least `min_size` elements.", call. = FALSE)
  } else if (length(keep_sets) < length(gene_sets)) {
    gene_sets <- gene_sets[keep_sets]
    set_sizes <- fsubset(set_sizes, keep_sets)
  }

  genes <- vec(gene_sets)
  unique_genes <- funique(genes)

  # Determine if any genes have an expected direction of change
  directional_sets <- anyv(grepl(";[ud]{1}", unique_genes, perl = TRUE), TRUE)

  if (directional_sets) {
    # Determine which genes are expected to be "down" and remove the direction
    # of change from the gene names
    gene_indices <- fmatch(genes, unique_genes)

    direction_down <- fsubset(
      .x = grepl(";d$", unique_genes, perl = TRUE),
      subset = gene_indices
    )

    unique_genes <- sub(";[ud]{1}$", "", unique_genes)

    genes <- fsubset(unique_genes, gene_indices)
  }

  # names(stats) is first because it was sorted lexicographically to deal with
  # ties in stats
  unique_genes <- intersect(names(stats), unique_genes)

  if (length(unique_genes) == 0L) {
    stop(
      "No elements of `gene_sets` are present in names(stats).",
      call. = FALSE
    )
  }

  # Indices of the genes in each set
  gene_indices <- fmatch(genes, unique_genes, nomatch = NA_integer_)

  # Unique integers for each set
  set_indices <- .C_rep_int(set_sizes, length(genes))

  unique_sets <- names(gene_sets)

  # Remove genes not in names(stats)
  if (anyNA(gene_indices)) {
    keep_genes <- whichNA(gene_indices, invert = TRUE)

    gene_indices <- fsubset(gene_indices, keep_genes)
    set_indices <- fsubset(set_indices, keep_genes)

    unique_sets <- fsubset(unique_sets, funique(set_indices))

    set_indices <- groupid(set_indices)
    class(set_indices) <- "integer"

    if (directional_sets) {
      direction_down <- fsubset(direction_down, keep_genes)
    }
  }

  n_sets <- length(unique_sets)

  # Split index vectors by direction of change
  if (directional_sets) {
    idx_down <- whichv(direction_down, TRUE)

    if (length(idx_down)) {
      gene_indices_down <- fsubset(gene_indices, idx_down)
      set_indices_down <- fsubset(set_indices, idx_down)

      gene_indices <- fsubset(gene_indices, -idx_down)
      set_indices <- fsubset(set_indices, -idx_down)
    } else {
      gene_indices_down <- set_indices_down <- integer(0L)
    }
  } else {
    gene_indices_down <- set_indices_down <- NULL
  }

  out <- list(
    "n_sets" = n_sets,
    "unique_genes" = unique_genes,
    "unique_sets" = unique_sets,
    "gene_indices" = gene_indices,
    "set_indices" = set_indices,
    "directional_sets" = directional_sets,
    "gene_indices_down" = gene_indices_down,
    "set_indices_down" = set_indices_down
  )

  return(out)
}


#' @title Fast, specialized `rep.int`
#'
#' @description Equivalent to `rep.int(seq_along(times), times)`, but several
#'   times faster when the output is large.
#'
#' @param times integer vector of group sizes.
#' @param length integer; the length of the result. Equal to `sum(times)`.
#'
#' @returns Integer vector.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.C_rep_int <- function(sizes, length) {
  .Call("_C_rep_int", sizes, length)
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


#' @title Remove gene sets that are too small or too large
#'
#' @param gene_indices integer vector; indices of the genes in each set,
#'   arranged in contiguous blocks by gene set.
#' @param extreme_set_indices integer vector; indices of gene sets that are too
#'   small or too large to test. Indices can range from 1 to `length(m)`.
#' @param m integer vector; the number of genes in each set, where
#'   `length(gene_indices) == sum(m)`.
#'
#' @returns The vector `gene_indices` with blocks of genes from the extreme gene
#'   sets removed.
#'
#' @details This function is not called when all or no gene sets are extreme. If
#'   all gene sets are extreme, an error will be thrown by `.calc_ES`.
#'
#'   The runtime of this function decreases as the number of extreme sets
#'   increases.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.C_remove_extreme_gene_sets <- function(gene_indices,
                                        extreme_set_indices,
                                        m) {
  .Call(
    "_C_remove_extreme_gene_sets",
    gene_indices,
    extreme_set_indices,
    m
  )
}


#' @title Calculate Enrichment Scores
#'
#' @param y_prime numeric vector of absolute gene-level values raised to the
#'   power of `alpha` for genes that are members of at least one gene set.
#' @param r_prime numeric vector of the ranks of the gene-level values for genes
#'   that are members of at least one gene set.
#' @param sum_ranks numeric; the sum of all ranks.
#' @param gene_indices integer vector; indices of the genes in all sets. Used to
#'   index vectors `y_prime` and `r_prime`.
#' @param m integer vector; the number of genes in each set. Used to select
#'   elements of `gene_indices`.
#' @param w integer vector; the number of genes that are not in each set.
#' @inheritParams hpgsea
#'
#' @returns Numeric vector of enrichment scores with the same length as `m`. If
#'   all elements of `y_prime` are 0 for a particular set, the ES for that set
#'   will be `NA`. If any `m` are less than `min_size`, the correspond ES will
#'   be 0.
#'
#' @author Tyler Sagendorf
#'
#' @noRd
.C_calc_ES <- function(y_prime,
                       r_prime,
                       sum_ranks,
                       gene_indices,
                       m,
                       w,
                       min_size) {
  .Call(
    "_C_calc_ES",
    y_prime,
    r_prime,
    sum_ranks,
    gene_indices,
    m,
    w,
    min_size
  )
}


#' @title Calculate Enrichment Scores and Other Information
#'
#' @description Calculate enrichment scores and information needed for
#'   permutation tests later on.
#'
#' @inheritParams hpgsea
#' @param n_genes integer; total number of genes. The length of `stats`.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom collapse fmatch fmax fsubset
#' @importFrom data.table frank
#'
#' @noRd
.calc_ES <- function(stats,
                     alpha = 1,
                     n_genes,
                     gene_sets,
                     min_size = 2L,
                     max_size = Inf) {
  max_size <- min(max_size, n_genes - 1L)

  storage.mode(min_size) <- storage.mode(max_size) <- "integer"

  index_list <- .gene_sets_to_indices(
    stats = stats,
    gene_sets = gene_sets,
    min_size = min_size
  )

  list2env(index_list, envir = environment())

  # Number of genes in each set, or the number of genes expected to be
  # up-regulated in each set
  m <- .C_group_sizes(set_indices, n_sets)

  if (directional_sets) {
    # Number of genes expected to be down-regulated in each set
    m_d <- .C_group_sizes(set_indices_down, n_sets)

    extreme_set_indices <- which(
      (m < min_size & m_d < min_size) | (m + m_d) > max_size
    )
  } else {
    m_d <- NULL

    extreme_set_indices <- which(
      m < min_size | m > max_size
    )
  }

  # Remove gene sets that are too small or too large
  if (length(extreme_set_indices) == n_sets) {
    stop(
      "All sets in `gene_sets` have fewer than `min_size` ",
      "or more than `max_size` genes in `stats`.",
      call. = FALSE
    )
  } else if (length(extreme_set_indices)) {
    unique_sets <- fsubset(unique_sets, -extreme_set_indices)

    gene_indices <- .C_remove_extreme_gene_sets(
      gene_indices = gene_indices,
      extreme_set_indices = extreme_set_indices,
      m = m
    )

    m <- fsubset(m, -extreme_set_indices)

    if (directional_sets) {
      gene_indices_down <- .C_remove_extreme_gene_sets(
        gene_indices = gene_indices_down,
        extreme_set_indices = extreme_set_indices,
        m = m_d
      )

      m_d <- fsubset(m_d, -extreme_set_indices)
    }
  }

  # The number of genes not in each set, or the number of genes that are not in
  # the set and expected to be up-regulated
  w <- n_genes - m

  if (directional_sets) {
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
  storage.mode(r) <- "double"

  # Subset to the genes that are elements of at least one gene set
  genes_in_sets <- fmatch(unique_genes, names(y))
  y_prime <- fsubset(y, genes_in_sets)
  r_prime <- fsubset(r, genes_in_sets)

  if (directional_sets) {
    # Calculate enrichment scores separately for the up-regulated and
    # down-regulated genes. Elements of m and m_d that are < min_size will be
    # replaced with 0.
    ES_u <- .C_calc_ES(
      y_prime = y_prime,
      r_prime = r_prime,
      sum_ranks = sum_ranks,
      gene_indices = gene_indices,
      m = m,
      w = w,
      min_size = min_size
    )

    ES_d <- .C_calc_ES(
      y_prime = y_prime,
      r_prime = r_prime,
      sum_ranks = sum_ranks,
      gene_indices = gene_indices_down,
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
      gene_indices = gene_indices,
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
    "w_d" = w_d,
    "directional_sets" = directional_sets
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
#' @author Pairing function devised by Matthew Szudzik. Code written by Tyler
#'   Sagendorf.
#'
#' @references Szudzik, M. (2006). An Elegant Pairing Function. Wolfram Science
#'   Conference. \url{http://szudzik.com/ElegantPairing.pdf}
#'
#' @noRd
.C_pair_szudzik <- function(x, y) {
  .Call("_C_pair_szudzik", x, y)
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
    unique_m <- funique(m)
    unique_w <- n_genes - unique_m

    map_unique_to_all_sets <- fmatch(m, unique_m)

    ES_end <- which(!duplicated(map_unique_to_all_sets, fromLast = TRUE))

    out <- list(
      "unique_m" = unique_m,
      "unique_w" = unique_w,
      "ES_end" = ES_end
    )
  } else {
    # Unique combinations of the number of up and down genes. Sorted by the
    # number of up-regulated genes and then by the number of down-regulated
    # genes.
    unique_pairs <- unique(cbind(m, m_d))

    unique_pairs_m_up <- unique_pairs[, 1L]
    unique_pairs_w_up <- n_genes - unique_pairs_m_up

    unique_pairs_m_down <- unique_pairs[, 2L]
    unique_pairs_w_down <- n_genes - unique_pairs_m_down

    # Unique number of up-regulated genes and the unique number of
    # down-regulated genes (separate, sorted vectors)
    unique_m_up <- funique(unique_pairs_m_up)
    unique_w_up <- n_genes - unique_m_up

    unique_m_down <- funique(unique_pairs_m_down, sort = TRUE)
    unique_w_down <- n_genes - unique_m_down

    map_pairs_to_all_sets <- fmatch(
      .C_pair_szudzik(m, m_d),
      .C_pair_szudzik(unique_pairs_m_up, unique_pairs_m_down)
    )

    map_unique_to_pairs_up <- fmatch(unique_pairs_m_up, unique_m_up)
    map_unique_to_pairs_down <- fmatch(unique_pairs_m_down, unique_m_down)

    ES_end <- which(!duplicated(map_pairs_to_all_sets, fromLast = TRUE))

    out <- list(
      "unique_m_up" = unique_m_up,
      "unique_w_up" = unique_w_up,
      "unique_m_down" = unique_m_down,
      "unique_w_down" = unique_w_down,
      "map_unique_to_pairs_up" = map_unique_to_pairs_up,
      "map_unique_to_pairs_down" = map_unique_to_pairs_down,
      "ES_end" = ES_end
    )
  }

  return(out)
}


#' @title Calculate Permutation Enrichment Scores
#'
#' @inheritParams hpgsea
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

  if (directional_sets) {
    tab[, set_size := m + m_d]
  } else {
    tab[, set_size := m]
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

  if (length(idx_not_NA) == 0L) {
    nperm <- 0L
  } else if (length(idx_not_NA) < length(ES)) {
    m <- fsubset(m, idx_not_NA)
    ES <- fsubset(ES, idx_not_NA)

    if (directional_sets) {
      m_d <- fsubset(m_d, idx_not_NA)
    }
  }

  # Vectors updated by permutation tests
  n_same_sign <- alloc(0L, length(ES))
  n_as_extreme <- alloc(0L, length(ES))
  sum_ES_perm <- alloc(0, length(ES))

  if (nperm > 0L) {
    size_list <- .unique_set_sizes(
      n_genes = n_genes,
      m = m,
      m_d = m_d
    )

    list2env(x = size_list, envir = environment())

    if (directional_sets) {
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
        unique_m_up,
        unique_w_up,
        unique_m_down,
        unique_w_down,
        map_unique_to_pairs_up,
        map_unique_to_pairs_down
      )
    } else {
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
        unique_m,
        unique_w
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
#' @inheritParams hpgsea
#'
#' @returns A `data.frame`.
#'
#' @author Tyler Sagendorf
#'
#' @importFrom collapse whichNA whichv
#' @importFrom data.table := setorderv setDF
#' @importFrom stats p.adjust
#'
#' @noRd
.calc_pvals <- function(tab,
                        nperm = 1e5L,
                        sort = TRUE,
                        alternative = "two.sided") {
  tab[, NES := ES / (sum_ES_perm / n_same_sign)]

  switch(
    EXPR = alternative,
    two.sided = {
      tab[
        whichNA(ES, invert = TRUE),
        p_value := (n_as_extreme + 1L) / (n_same_sign + 1L)
      ]
    },
    less = {
      tab[
        whichNA(ES, invert = TRUE),
        p_value := ifelse(
          ES < 0,
          n_as_extreme + 1L,
          nperm - n_as_extreme + 1L
        ) / (nperm + 1L)
      ]
    },
    greater = {
      tab[
        whichNA(ES, invert = TRUE),
        p_value := ifelse(
          ES >= 0,
          n_as_extreme + 1L,
          nperm - n_as_extreme + 1L
        ) / (nperm + 1L)
      ]
    }
  )

  tab[
    whichv(n_same_sign, 0L), # only happens when nperm is very small
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
