test_that(".C_remove_extreme_sets is correct", {
  set_sizes <- sample.int(50L, size = 1000L, replace = TRUE)

  sum_sizes <- sum(set_sizes)

  gene_indices <- seq_len(sum_sizes)
  set_indices <- .C_rep_int(set_sizes, length = sum_sizes)

  extreme_set_indices <- which(
    set_sizes < 5L | set_sizes > 40L
  )

  res <- .C_remove_extreme_gene_sets(
    gene_indices = gene_indices,
    extreme_set_indices = extreme_set_indices,
    m = set_sizes
  )

  expected <- fsubset(
    .x = gene_indices,
    # collapse::`%!iin%`
    subset = whichNA(fmatch(set_indices, extreme_set_indices))
  )

  expect_identical(
    res,
    expected
  )
})
