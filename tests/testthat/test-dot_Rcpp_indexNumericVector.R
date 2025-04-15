test_that(".Rcpp_indexNumericVector works correctly", {
   set.seed(0)

   x <- rnorm(10L)
   idx <- dqsample.int(n = 10L, size = 200L, replace = TRUE)

   expected <- x[idx]

   actual <- .Rcpp_indexNumericVector(x, idx)

   expect_identical(object = actual,
                    expected = expected)
})
