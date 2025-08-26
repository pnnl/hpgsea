test_that("read_gmt creates a named list of sets from a GMT file", {
  expect_no_error({
    gene_sets <- read_gmt(file = gmt_file)
  })

  # Should be TRUE
  res <- is.list(gene_sets) &&
    !is.null(names(gene_sets)) &&
    length(gene_sets) == 50L

  expect_true(res)
})
