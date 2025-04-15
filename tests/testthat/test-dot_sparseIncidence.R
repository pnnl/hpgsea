
test_that("x is a named list of character vectors", {
   err1 <- capture_error(
      .sparseIncidence(x = letters)
   )$message

   err2 <- capture_error(
      .sparseIncidence(x = list("Set1" = 1:3))
   )$message

   expected_err <- "`x` must be a named list of character vectors."

   expect_identical(
      object = err1,
      expected = expected_err
   )

   expect_identical(
      object = err2,
      expected = expected_err
   )
})


test_that("the sets overlap with the background", {
   background <- LETTERS

   err <- capture_error(
      .sparseIncidence(x = list("Set1" = letters),
                       background = LETTERS)
   )$message

   expect_identical(
      object = err,
      expected = "No elements of `gene_sets` are present in `background`."
   )
})



test_that("The non-directional incidence matrix is correct", {
   # Create gene sets
   gene_sets <- list("Set2" = c("a", "b", "d", "h", "d"),
                     "Set1" = c("d", "f", "n"))

   # Shuffle the elements. The rows of the incidence matrix should be in
   # alphabetical order, so this should not affect the results.
   gene_sets <- lapply(gene_sets, sample)

   # This is what the incidence matrix should look like
   A_true <- Matrix::sparseMatrix(i = c(1L, 2L, 3L, 3L, 4L, 5L),
                                  j = c(1L, 1L, 1L, 2L, 2L, 1L),
                                  x = rep(1, 6L),
                                  dims = c(5L, 2L),
                                  dimnames = list(c("a", "b", "d", "f", "h"),
                                                  c("Set2", "Set1")))

   # The background has duplicates, which should not affect the results
   A_list <- .sparseIncidence(x = gene_sets,
                              background = letters[c(1:10, 10:5)])

   A <- A_list[["A"]]
   A.d <- A_list[["A.d"]]

   # Checks
   expect_identical(object = A, expected = A_true)

   expect_null(A.d) # not a directional database
})


test_that("The directional incidence matrix is correct", {
   # Create gene sets (include duplicate elements)
   gene_sets <- list("Set3" = c("a;u", "a;u", "b;d", "c;u"),
                     "Set1" = c("r;d", "b;u", "c;d"),
                     "Set2" = c("z;d", "z;d"))

   gene_sets <- lapply(gene_sets, sample) # shuffle elements

   ## Expected incidence matrices
   A_true <- Matrix::sparseMatrix(
      i = c(1L, 2L, 3L),
      j = c(1L, 2L, 1L),
      x = rep(1, 3L),
      dims = c(5L, 3L),
      dimnames = list(c("a", "b", "c", "r", "z"),
                      c("Set3", "Set1", "Set2"))
   )

   A_d_true <- Matrix::sparseMatrix(
      i = c(2L, 3L, 4L, 5L),
      j = c(1L, 2L, 2L, 3L),
      x = rep(1, 4L),
      dims = c(5L, 3L),
      dimnames = list(c("a", "b", "c", "r", "z"),
                      c("Set3", "Set1", "Set2"))
   )

   # Includes duplicates
   background <- c("a", "b", "c", "c", "r", "z", "t")

   # List of incidence matrices
   A_list <- .sparseIncidence(x = gene_sets,
                              background = background)

   A <- A_list[["A"]] # genes expected to be "up" (;u)
   A_d <- A_list[["A_d"]] # genes expected to be "down" (;d)

   # Checks
   expect_identical(object = A, expected = A_true)

   expect_identical(object = A_d, expected = A_d_true)

   ## Elements can not be both "up" and "down" in the same set
   gene_sets[["Set3"]] <- c(gene_sets[["Set3"]], "a;d") # a;u already exists

   err <- capture_error(
      fast.ssgsea:::.sparseIncidence(x = gene_sets,
                                     background = background),
   )$message

   expect_identical(
      object = err,
      expected = paste0("Elements can not be both up (suffix \";u\") and ",
                        "down (suffix \";d\") in the same set.")
   )

})
