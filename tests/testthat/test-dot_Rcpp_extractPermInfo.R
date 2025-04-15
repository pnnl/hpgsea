test_that("All values are 0 when signs mismatch", {
   x <- c(0, 2, 4.5, 6.7, 12)
   y <- c(-10, -8, -6, -2)

   expected <- list(
      "n_same_sign_b" = rep(0L, 5L),
      "n_as_extreme_b" = rep(0L, 5L),
      "sum_ES_perm_b" = rep(0, 5L)
   )
   actual <- .Rcpp_extractPermInfo(x, y)

   expect_identical(expected, actual)

   x <- c(-14, -10.8, -6, -4, -0.1)
   y <- c(0, 3, 5, 8.2, 1, 4, 2, 2)

   actual <- .Rcpp_extractPermInfo(x, y)

   expect_identical(expected, actual)
})


test_that("Positive results are correct", {
   x <- c(1, 2.2, 4.5, 7, 12)
   y <- c(4.5, -11.2, 2.8, -13, -2, 4.7, 8, 10, 6.1)

   expected <- list(
      "n_same_sign_b" = rep(6L, 5L),
      "n_as_extreme_b" = c(6L, 6L, 5L, 2L, 0L),
      "sum_ES_perm_b" = rep(sum(y[y >= 0]), 5L)
   )

   actual <- .Rcpp_extractPermInfo(x, y)

   expect_identical(expected, actual)
})


test_that("Negative results are correct", {
   x <- rev(-c(1, 2.2, 4.5, 7, 12))
   y <- c(4.5, -11.2, 2.8, -13, -2, 4.7, 8, 10, 6.1)

   expected <- list(
      "n_same_sign_b" = rep(3L, 5L),
      "n_as_extreme_b" = c(1L, 2L, 2L, 2L, 3L),
      "sum_ES_perm_b" = rep(-sum(y[y < 0]), 5L)
   )

   actual <- .Rcpp_extractPermInfo(x, y)

   expect_identical(expected, actual)
})


test_that("Mixed sign results are correct", {
   x <- c(-12, -8, -6, -4, 0, 1, 5, 20)
   y <- sample(c(-10, -8, -7, -4, -3, 0, 1, 4, 5, 7, 8, 15))

   expected <- list(
      "n_same_sign_b" = rep(c(5L, 7L), each = 4L),
      "n_as_extreme_b" = c(0L, 2L, 3L, 4L, 7L, 6L, 4L, 0L),
      "sum_ES_perm_b" = rep(c(32, 40), each = 4L)
   )

   actual <- .Rcpp_extractPermInfo(x, y)

   expect_identical(expected, actual)
})
