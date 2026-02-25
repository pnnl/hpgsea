# These tests are a mess. Some of them can be moved to separate scripts for
# their corresponding internal function.

test_that("nperm <= 2 billion", {
  err <- capture_error(
    fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = gene_sets,
      nperm = 2e9L + 1L
    )
  )$message

  expect_identical(
    object = err,
    expected = "`nperm` must be an integer between 0 and 2 billion."
  )
})


test_that("alpha is finite and non-negative", {
  err <- capture_error(
    fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = gene_sets,
      alpha = Inf
    )
  )$message

  expect_identical(
    object = err,
    expected = "`alpha` must be a single non-negative real number."
  )
})


test_that("`sort` must be logical", {
  expected_error <- "`sort` must be TRUE or FALSE."

  err1 <- capture_error(
    fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = gene_sets,
      sort = list(logical(1L))
    )
  )$message

  err2 <- capture_error(
    fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = gene_sets,
      sort = 1
    )
  )$message

  expect_identical(err1, expected_error)
  expect_identical(err2, expected_error)
})


test_that("`stats` is a numeric vector with unique names", {
  err1 <- capture_error(
    fast_ssgsea(
      stats = numeric(1L),
      gene_sets = list("Set1" = letters)
    )
  )$message

  err2 <- capture_error(
    fast_ssgsea(
      stats = structure(c(1, 2), names = c("a", "a")),
      gene_sets = list("Set1" = letters)
    )
  )$message

  expect_identical(err1, err2)

  expect_identical(
    object = err1,
    expected = "`stats` must be a numeric vector with unique names."
  )
})


test_that("`stats` has at least 3 nonmissing values", {
  err1 <- capture_error(
    fast_ssgsea(
      stats = stats_mat[1:2, 1L],
      gene_sets = gene_sets
    )
  )$message

  expect_identical(
    object = err1,
    expected = "`stats` must have at least 3 nonmissing values."
  )
})


test_that("`min_size` is smaller than the number of nonmissing values", {
  err1 <- capture_error(
    fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = gene_sets,
      nperm = 0L,
      min_size = nrow(stats_mat)
    )
  )$message

  expect_identical(
    object = err1,
    expected = paste0(
      "`min_size` must be >= 2 and less than the number of ",
      "nonmissing values in `stats`."
    )
  )
})


test_that("extreme sets are removed", {
  # Remove extreme sets, unless all sets will be removed
  res <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = c(gene_sets, list("Set999" = rownames(stats_mat))),
    nperm = 0L
  )

  expect_false("Set999" %in% res$GeneSet)

  # If all sets are extreme, throw an error
  err <- capture_error(
    fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = list(
        "Set1" = rownames(stats_mat),
        "Set2" = rownames(stats_mat)[1:2]
      ),
      nperm = 0L,
      min_size = 3L
    )
  )$message

  expect_identical(
    object = err,
    expected = paste0(
      "All sets in `gene_sets` have fewer than `min_size` or ",
      "more than `max_size` genes."
    )
  )
})


test_that("extreme directional sets are removed", {
  # Remove extreme sets, unless all sets will be removed
  res <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = c(
      gene_sets,
      list(
        "Set999" = paste0(rownames(stats_mat), ";d"),
        "Set998" = paste0(rownames(stats_mat)[1L], ";u")
      )
    ),
    nperm = 0L
  )

  expect_true(
    length(intersect(c("Set999", "Set998"), res$GeneSet)) == 0L
  )

  # If all sets are extreme, throw an error
  err <- capture_error(
    fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = list(
        "Set1" = paste0(rownames(stats_mat), ";u"),
        "Set2" = paste0(rownames(stats_mat)[1:2], ";d")
      ),
      nperm = 0L,
      min_size = 3L
    )
  )$message

  expect_identical(
    object = err,
    expected = paste0(
      "All sets in `gene_sets` have fewer than `min_size` or ",
      "more than `max_size` genes."
    )
  )
})


test_that("The columns of the results have the correct type and position", {
  res <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    alpha = 0,
    nperm = 100L,
    sort = FALSE
  )

  expected_cols <- c(
    "sample", "set", "set_size",
    "ES", "NES", "n_same_sign",
    "n_as_extreme", "p_value", "adj_p_value"
  )

  expect_identical(
    object = colnames(res),
    expected = expected_cols
  )

  expected_col_types <- c(
    "factor", "character", "integer",
    "numeric", "numeric", "integer",
    "integer", "numeric", "numeric"
  )
  names(expected_col_types) <- expected_cols

  expect_identical(
    object = vapply(res, class, character(1L)),
    expected = expected_col_types
  )

  ## Directional database
  set.seed(0)
  gene_set1_up <- paste0(
    sample(rownames(stats_mat)[1:150], size = 60L),
    ";u"
  )
  gene_set1_down <- paste0(
    sample(rownames(stats_mat)[nrow(stats_mat) - 1:150], size = 40L),
    ";d"
  )

  gene_sets_dir <- list("Set1" = c(gene_set1_up, gene_set1_down))

  res <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets_dir,
    nperm = 10L
  )

  col_idx <- match(c("ES_u", "ES_d", "ES"), colnames(res))

  expect_identical(
    object = col_idx,
    expected = 4:6
  )
})


test_that("The ES are correct when there are ties", {
  ## alpha = 0
  expected_ES_0 <- lapply(gene_sets, function(set_i) {
    calculate_ES(stats_mat = stats_mat, gene_set = set_i, alpha = 0)
  })
  expected_ES_0 <- data.table::rbindlist(expected_ES_0, idcol = "set")
  setorderv(expected_ES_0, cols = c("sample", "set"))

  true_ES_0 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    alpha = 0,
    nperm = 0L,
    sort = FALSE
  )

  expect_equal(
    object = true_ES_0$ES,
    expected = expected_ES_0$ES
  )


  ## alpha = 1
  expected_ES_1 <- lapply(gene_sets, function(set_i) {
    calculate_ES(stats_mat = stats_mat, gene_set = set_i, alpha = 1)
  })
  expected_ES_1 <- data.table::rbindlist(expected_ES_1, idcol = "set")
  setorderv(expected_ES_1, cols = c("sample", "set"))

  true_ES_1 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    alpha = 1,
    nperm = 0L,
    sort = FALSE
  )

  expect_equal(
    object = true_ES_1$ES,
    expected = expected_ES_1$ES
  )
})


test_that("the ES are correct for directional sets", {
  # Directional sets
  set.seed(0)

  # First gene set. All genes have extreme values in the tail of their expected
  # direction of change, indicating really good agreement. The ES and NES
  # should be very positive.
  gene_set1_up <- paste0(
    sample(rownames(stats_mat)[1:150], size = 60L),
    ";u"
  )
  gene_set1_down <- paste0(
    sample(rownames(stats_mat)[nrow(stats_mat) - 1:150], size = 40L),
    ";d"
  )

  # Genes in set 2 are only "down", but the values are positive
  gene_set2_down <- paste0(
    sample(rownames(stats_mat)[1:300], size = 30L),
    ";d"
  )

  gene_sets_dir <- list(
    "Set1" = c(gene_set1_up, gene_set1_down),
    "Set2" = gene_set2_down,
    "Set3" = gene_set1_up
  )

  ## alpha = 1 ----
  res1 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets_dir,
    alpha = 1,
    nperm = 1e4L,
    sort = FALSE,
    seed = 0
  )

  # Need to calculate ES in a piece-wise fashion for directional sets
  set1_res1 <- calculate_ES(
    stats_mat = stats_mat,
    gene_set = sub(";u$", "", gene_set1_up),
    alpha = 1
  )
  colnames(set1_res1) <- c("sample", "ES_u")

  set1_res2 <- calculate_ES(
    stats_mat = stats_mat,
    gene_set = sub(";d$", "", gene_set1_down),
    alpha = 1
  )
  colnames(set1_res2) <- c("sample", "ES_d")

  set2_res <- calculate_ES(
    stats_mat = stats_mat,
    gene_set = sub(
      ";[ud]$", "", gene_sets_dir[["Set2"]]
    ),
    alpha = 1
  )
  set2_res[, `:=`(set = "Set2", ES_u = 0, ES_d = ES)][, ES := -ES]
  set2_res <- set2_res[, c("sample", "set", "ES_u", "ES_d", "ES")]

  set3_res <- calculate_ES(
    stats_mat = stats_mat,
    gene_set = sub(
      ";[ud]$", "", gene_sets_dir[["Set3"]]
    ),
    alpha = 1
  )
  set3_res[, `:=`(set = "Set3", ES_d = 0, ES_u = ES)]
  set3_res <- set3_res[, c("sample", "set", "ES_u", "ES_d", "ES")]

  set1_res <- merge(x = set1_res1, y = set1_res2, by = "sample")
  set1_res[, `:=`(set = "Set1", ES = ES_u - ES_d)]
  data.table::setcolorder(set1_res, neworder = "set", after = "sample")

  true_res1 <- rbind(set1_res, set2_res, set3_res)
  setorderv(true_res1, cols = c("sample", "set"))
  true_res1 <- as.data.frame(true_res1)

  expect_equal(
    object = res1[, colnames(true_res1)],
    expected = true_res1
  )

  # The "up" enrichment scores for Set2 should be 0, since no genes in the set
  # had the ";u" suffix.
  expect_true(
    all(res1$ES_u[res1$set == "Set2"] == 0)
  )

  # The "down" enrichment scores for Set3 should be 0, since no genes in the
  # set had the ";d" suffix.
  expect_true(
    all(res1$ES_d[res1$set == "Set3"] == 0)
  )

  # No NES should be NA
  expect_true(
    all(!is.na(res1$NES))
  )

  # The gene set is the same as gene_set2_down (without the expected direction
  # of change), but the sign of the ES and NES will have opposite signs to the
  # res1 results.
  res2 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = list("Set1" = sub(";d", "", gene_set2_down)),
    alpha = 1,
    nperm = 1e4L,
    sort = FALSE,
    seed = 0
  )

  expect_equal(
    res1$ES[res1$set == "Set2"],
    -1 * res2$ES
  )

  # Slight differences due to using floats instead of doubles to calculate
  # permutation enrichment scores in Rcpp_calcESPermCore()
  expect_equal(
    signif(res1$NES[res1$set == "Set2"], digits = 5L),
    c(-6.15990, -0.16928, 0.77529, 0.09746)
  )


  ## alpha = 0 ----
  res1 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets_dir,
    alpha = 0,
    nperm = 1e4L,
    sort = FALSE,
    seed = 0
  )

  # The "up" enrichment scores for Set2 should be 0, since no genes in the set
  # had the ";u" suffix.
  expect_true(
    all(res1$ES_u[res1$set == "Set2"] == 0)
  )

  # No NES should be NA
  expect_true(
    all(!is.na(res1$NES))
  )

  # The gene set is the same as gene_set2_down (without the expected direction
  # of change), but the sign of the ES and NES will have opposite signs to the
  # res1 results.
  res2 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = list("Set1" = sub(";d", "", gene_set2_down)),
    alpha = 0,
    nperm = 1e4L,
    sort = FALSE,
    seed = 0
  )

  expect_equal(
    res1$ES[res1$set == "Set2"],
    -1 * res2$ES
  )

  expect_equal(
    signif(res1$NES[res1$set == "Set2"], digits = 5L),
    c(-10.41400, 0.064937, 0.374570, -0.050986)
  )
})


test_that("genes in directional gene sets can be the same direction", {
  # All down
  set_down <- paste0(rownames(stats_mat)[seq_len(30L)], ";d")

  expect_no_error(
    res1 <- fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = list("set_down" = set_down),
      nperm = 500L
    )
  )

  expect_true(
    all(!is.na(res1$NES))
  )

  # All up
  set_up <- paste0(rownames(stats_mat)[seq_len(30L)], ";u")

  expect_no_error(
    res2 <- fast_ssgsea_multicol(
      stats_mat = stats_mat,
      gene_sets = list("set_up" = set_up),
      nperm = 500L
    )
  )

  expect_true(
    all(!is.na(res2$NES))
  )

})


test_that("NES are mostly within [-4, +4]", {
  set_sizes <- rep.int(5:100, 20L)

  gene_sets <- lapply(seq_along(set_sizes), function(i) {
    set.seed(i)

    sample(rownames(stats_mat), size = set_sizes[i])
  })
  names(gene_sets) <- paste0("set.", seq_along(gene_sets))

  res1 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    alpha = 0,
    nperm = 500L,
    sort = FALSE,
    seed = 0L
  )

  expect_true(
    mean(res1$NES <= 4 & res1$NES >= -4) >= 0.995
  )

  res2 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    alpha = 1,
    nperm = 500L,
    sort = FALSE,
    seed = 0L
  )

  expect_true(
    mean(res2$NES <= 4 & res2$NES >= -4) >= 0.995
  )
})


test_that("results are sorted correctly", {
  res1 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    alpha = 0,
    nperm = 500L,
    sort = FALSE
  )

  # Gene sets should appear in order (no sorting)
  expect_identical(
    object = res1$set,
    expected = rep(names(gene_sets), times = ncol(stats_mat))
  )
})


test_that("n_same_sign >= n_as_extreme", {
  res <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    nperm = 500L,
    sort = FALSE
  )

  expect_true(
    all(res$n_same_sign >= res$n_as_extreme)
  )
})


test_that("p-values are between 0 and 1", {
  res1 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    alpha = 0,
    nperm = 500L,
    sort = FALSE
  )

  expect_true(
    all(res1$p_value <= 1) && all(res1$p_value > 0)
  )

  res2 <- fast_ssgsea_multicol(
    stats_mat = stats_mat,
    gene_sets = gene_sets,
    alpha = 1,
    nperm = 500L,
    sort = FALSE
  )

  expect_true(
    all(res2$p_value <= 1) && all(res2$p_value > 0)
  )
})
