test_that(".Rcpp_matmult_dense works correctly", {
   set.seed(0L)

   x <- rnorm(3 * 3)
   y <- rnorm(3 * 7)

   dim(x) <- c(3L, 3L)
   dim(y) <- c(3L, 7L)

   expected <- x %*% y

   actual <- .Rcpp_matmult_dense(x, y)

   expect_identical(actual, expected)
})
