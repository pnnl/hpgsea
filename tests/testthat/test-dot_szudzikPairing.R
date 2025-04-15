test_that("Szudzik pairing function is correctly specified", {
  x <- seq_len(5L)
  y <- 5L - x + 1L

  out <- .szudzikPairing(x, y)

  expect_identical(object = out,
                   expected = c(26L, 18L, 15L, 22L, 31L))
})
