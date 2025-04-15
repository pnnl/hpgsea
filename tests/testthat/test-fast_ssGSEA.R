# These tests are a mess. Some of them can be moved to separate scripts for
# their corresponding internal function.

test_that("nperm <= 1 million", {
   err <- capture_error(
      fast_ssgsea(X = X, gene_sets = gene_sets, nperm = 1e6L + 1L)
   )$message

   expect_identical(
      object = err,
      expected = "`nperm` must be a whole number between 0 and 1 million."
   )
})


test_that("alpha is finite and non-negative", {
   err <- capture_error(
      fast_ssgsea(X = X, gene_sets = gene_sets, alpha = Inf)
   )$message

   expect_identical(
      object = err,
      expected = "`alpha` must be a single non-negative real number."
   )
})


test_that("`batch_size` is between 1 and `nperm`", {
   expected_error <- "`batch_size` must be a whole number between 1 and `nperm`."

   err1 <- capture_error(
      fast_ssgsea(X = X, gene_sets = gene_sets, batch_size = 0L)
   )$message

   err2 <- capture_error(
      fast_ssgsea(X = X, gene_sets = gene_sets, batch_size = 0.5)
   )$message

   err3 <- capture_error(
      fast_ssgsea(X = X, gene_sets = gene_sets, batch_size = NA_real_)
   )$message

   expect_identical(err1, expected_error)
   expect_identical(err2, expected_error)
   expect_identical(err3, expected_error)
})


test_that("`sort` must be logical", {
   expected_error <- "`sort` must be TRUE or FALSE."

   err1 <- capture_error(
      fast_ssgsea(X = X,
                  gene_sets = gene_sets,
                  sort = list(logical(1L)))
   )$message

   err2 <- capture_error(
      fast_ssgsea(X = X,
                  gene_sets = gene_sets,
                  sort = 1)
   )$message

   expect_identical(err1, expected_error)
   expect_identical(err2, expected_error)
})


test_that("`adjust_globally` must be logical", {
   expected_error <- "`adjust_globally` must be TRUE or FALSE."

   err1 <- capture_error(
      fast_ssgsea(X = X,
                  gene_sets = gene_sets,
                  adjust_globally = list(logical(1L)))
   )$message

   err2 <- capture_error(
      fast_ssgsea(X = X,
                  gene_sets = gene_sets,
                  adjust_globally = 1)
   )$message

   expect_identical(err1, expected_error)
   expect_identical(err2, expected_error)
})


test_that("`X` is a numeric matrix with row and column names", {
   err1 <- capture_error(
      fast_ssgsea(X = numeric(1L),
                  gene_sets = list("Set1" = letters))
   )$message

   err2 <- capture_error(
      fast_ssgsea(X = matrix(rnorm(4), nrow = 2L, ncol = 2L),
                  gene_sets = list("Set1" = letters))
   )$message

   err3 <- capture_error(
      fast_ssgsea(X = matrix(c(TRUE, FALSE), nrow = 2L, ncol = 1L),
                  gene_sets = list("Set1" = letters))
   )$message

   expect_true(length(unique(c(err1, err2, err3))) == 1L)

   expect_identical(
      object = err1,
      expected = "`X` must be a numeric matrix with row and column names."
   )
})


test_that("`X` has at least 3 rows", {
   err1 <- capture_error(
      fast_ssgsea(X = X[1:2, ], gene_sets = gene_sets)
   )$message

   expect_identical(
      object = err1,
      expected = "Matrix `X` must have at least 3 rows."
   )
})


test_that("`X` has enough nonmissing values in each sample", {
   X2 <- X
   X2[3:nrow(X), 2L] <- NA

   err <- capture_error(
      fast_ssgsea(X = X2,
                  gene_sets = gene_sets,
                  nperm = 0L)
   )$message

   expect_identical(
      object = err,
      expected = c(
         "Matrix `X` must have at least 3 nonmissing values in each column."
      )
   )
})


test_that("`min_size` is smaller than the number of nonmissing values", {
   err <- capture_error(
      fast_ssgsea(X = X,
                  gene_sets = gene_sets,
                  nperm = 0L,
                  min_size = nrow(X))
   )$message

   expect_identical(
      object = err,
      expected = "`min_size` must be >= 2 and < nrow(X)."
   )
})


test_that("extreme sets are removed", {
   # Remove extreme sets, unless all sets will be removed
   res <- fast_ssgsea(X = X,
                      gene_sets = c(gene_sets, list("Set999" = rownames(X))),
                      nperm = 0L)

   expect_false("Set999" %in% res$GeneSet)

   # If all sets are extreme, throw an error
   err <- capture_error(
      fast_ssgsea(X = X,
                  gene_sets = list(
                     "Set1" = rownames(X),
                     "Set2" = rownames(X)[1:2]
                  ),
                  nperm = 0L,
                  min_size = 3L)
   )$message

   expect_identical(
      object = err,
      expected = paste0(
         "All sets in `gene_sets` contain fewer than `min_size` genes with ",
         "nonmissing values or consist of all genes with nonmissing values ",
         "in at least one sample."
      )
   )
})


test_that("extreme directional sets are removed", {
   # Remove extreme sets, unless all sets will be removed
   res <- fast_ssgsea(X = X,
                      gene_sets = c(gene_sets,
                                    list("Set999" = paste0(rownames(X), ";d"),
                                         "Set998" = paste0(rownames(X)[1L], ";u"))),
                      nperm = 0L)

   expect_true(
      length(intersect(c("Set999", "Set998"), res$GeneSet)) == 0L
   )

   # If all sets are extreme, throw an error
   err <- capture_error(
      fast_ssgsea(X = X,
                  gene_sets = list(
                     "Set1" = paste0(rownames(X), ";u"),
                     "Set2" = paste0(rownames(X)[1:2], ";d")
                  ),
                  nperm = 0L,
                  min_size = 3L)
   )$message

   expect_identical(
      object = err,
      expected = paste0(
         "All sets in `gene_sets` contain fewer than `min_size` genes with ",
         "nonmissing values or consist of all genes with nonmissing values ",
         "in at least one sample."
      )
   )
})


test_that("The columns of the results have the correct type and position", {
   res <- fast_ssgsea(X = X,
                      gene_sets = gene_sets,
                      alpha = 0,
                      nperm = 100L,
                      sort = FALSE)

   expected_cols <- c("sample", "set", "set_size",
                      "ES", "NES", "n_same_sign",
                      "n_as_extreme", "p_value", "adj_p_value")

   expect_identical(object = colnames(res),
                    expected = expected_cols)

   expected_col_types <- c("factor", "character", "integer",
                           "numeric", "numeric", "integer",
                           "integer", "numeric", "numeric")
   names(expected_col_types) <- expected_cols

   expect_identical(object = vapply(res, class, character(1L)),
                    expected = expected_col_types)

   ## Directional database
   set.seed(0)
   gene_set1_up <- paste0(sample(rownames(X)[1:150], size = 60L),
                          ";u")
   gene_set1_down <- paste0(sample(rownames(X)[nrow(X) - 1:150], size = 40L),
                            ";d")

   gene_sets_dir <- list("Set1" = c(gene_set1_up, gene_set1_down))

   res <- fast_ssgsea(X = X,
                      gene_sets = gene_sets_dir,
                      nperm = 10L)

   col_idx <- match(c("ES_u", "ES_d", "ES"), colnames(res))

   expect_identical(object = col_idx,
                    expected = 4:6)
})


test_that("The ES are correct when there are ties", {
   ## alpha = 0
   expected_ES_0 <- lapply(gene_sets, function(set_i) {
      calculate_ES(X = X, gene_set = set_i, alpha = 0)
   })
   expected_ES_0 <- rbindlist(expected_ES_0, idcol = "set")
   setorderv(expected_ES_0, cols = c("sample", "set"))

   true_ES_0 <- fast_ssgsea(X = X,
                            gene_sets = gene_sets,
                            alpha = 0,
                            nperm = 0L,
                            sort = FALSE)

   expect_equal(object = true_ES_0$ES,
                expected = expected_ES_0$ES)


   ## alpha = 1
   expected_ES_1 <- lapply(gene_sets, function(set_i) {
      calculate_ES(X = X, gene_set = set_i, alpha = 1)
   })
   expected_ES_1 <- rbindlist(expected_ES_1, idcol = "set")
   setorderv(expected_ES_1, cols = c("sample", "set"))

   true_ES_1 <- fast_ssgsea(X = X,
                            gene_sets = gene_sets,
                            alpha = 1,
                            nperm = 0L,
                            sort = FALSE)

   expect_equal(object = true_ES_1$ES,
                expected = expected_ES_1$ES)
})


test_that("the ES are correct for directional sets", {
   # Directional sets
   set.seed(0)

   # First gene set. All genes have extreme values in the tail of their expected
   # direction of change, indicating really good agreement. The ES and NES
   # should be very positive.
   gene_set1_up <- paste0(sample(rownames(X)[1:150], size = 60L),
                          ";u")
   gene_set1_down <- paste0(sample(rownames(X)[nrow(X) - 1:150], size = 40L),
                            ";d")

   # Genes in set 2 are only "down", but the values are positive
   gene_set2_down <- paste0(sample(rownames(X)[1:300], size = 30L),
                            ";d")

   gene_sets_dir <- list("Set1" = c(gene_set1_up, gene_set1_down),
                         "Set2" = gene_set2_down,
                         "Set3" = gene_set1_up)

   ## alpha = 1 ----
   res1 <- fast_ssgsea(X = X,
                       gene_sets = gene_sets_dir,
                       alpha = 1,
                       nperm = 1000L,
                       sort = FALSE,
                       seed = 0)

   # Need to calculate ES in a piece-wise fashion for directional sets
   set1_res1 <- calculate_ES(X = X,
                             gene_set = sub(";u$", "", gene_set1_up),
                             alpha = 1)
   colnames(set1_res1) <- c("sample", "ES_u")

   set1_res2 <- calculate_ES(X = X,
                             gene_set = sub(";d$", "", gene_set1_down),
                             alpha = 1)
   colnames(set1_res2) <- c("sample", "ES_d")

   set2_res <- calculate_ES(X = X,
                            gene_set = sub(
                               ";[ud]$", "", gene_sets_dir[["Set2"]]
                            ),
                            alpha = 1)
   set2_res[, `:=`(set = "Set2", ES_u = 0, ES_d = ES)][, ES := -ES]
   set2_res <- set2_res[, c("sample", "set", "ES_u", "ES_d", "ES")]

   set3_res <- calculate_ES(X = X,
                            gene_set = sub(
                               ";[ud]$", "", gene_sets_dir[["Set3"]]
                            ),
                            alpha = 1)
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
   res2 <- fast_ssgsea(X = X,
                       gene_sets = list("Set1" = sub(";d", "", gene_set2_down)),
                       alpha = 1,
                       nperm = 1000L,
                       sort = FALSE,
                       seed = 0)

   expect_equal(
      res1$ES[res1$set == "Set2"],
      -1 * res2$ES
   )

   # Slight differences due to using floats instead of doubles to calculate
   # permutation enrichment scores in Rcpp_calcESPermCore()
   expect_equal(
      signif(res1$NES[res1$set == "Set2"], digits = 6L),
      signif(-1 * res2$NES, digits = 6L)
   )


   ## alpha = 0 ----
   res1 <- fast_ssgsea(X = X,
                       gene_sets = gene_sets_dir,
                       alpha = 0,
                       nperm = 1000L,
                       sort = FALSE,
                       seed = 0)

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
   res2 <- fast_ssgsea(X = X,
                       gene_sets = list("Set1" = sub(";d", "", gene_set2_down)),
                       alpha = 0,
                       nperm = 1000L,
                       sort = FALSE,
                       seed = 0)

   expect_equal(
      res1$ES[res1$set == "Set2"],
      -1 * res2$ES
   )

   expect_equal(
      res1$NES[res1$set == "Set2"],
      -1 * res2$NES
   )
})


test_that("results are sorted correctly", {
   res1 <- fast_ssgsea(X = X,
                       gene_sets = gene_sets,
                       alpha = 0,
                       nperm = 500L,
                       sort = FALSE)

   # Samples should appear one after the other
   expected_samples <- paste0("sample", rep(seq_len(ncol(X)), each = 2L))
   expected_samples <- as.factor(expected_samples)

   expect_identical(object = res1$sample,
                    expected = expected_samples)

   # Gene sets should appear in order (no sorting)
   expect_identical(object = res1$set,
                    expected = rep(names(gene_sets), times = ncol(X)))
})


test_that("p-values are adjusted separately by sample", {
   ## Adjust p-values separately by sample
   res1 <- fast_ssgsea(X = X,
                       gene_sets = gene_sets,
                       alpha = 0,
                       nperm = 500L,
                       sort = FALSE,
                       adjust_globally = FALSE)

   expect_identical(
      object = res1$adj_p_value[seq_along(gene_sets)],
      expected = p.adjust(res1$adj_p_value[seq_along(gene_sets)],
                          method = "BH")
   )

   ## Adjust p-values separately by sample
   res2 <- fast_ssgsea(X = X,
                       gene_sets = gene_sets,
                       nperm = 500L,
                       sort = FALSE,
                       adjust_globally = TRUE)

   expect_identical(
      object = res2$adj_p_value,
      expected = p.adjust(res2$p_value, method = "BH")
   )
})
