# Write a named list of gene sets to a GMT file. Faster than cmapR version.
writeGMT <- function(x, path) {

   if (!grepl("\\.gmt$", path))
      stop("`path` is not a path to a GMT file.")

   x <- vapply(x, function(xi) {
      paste(xi, collapse = "\t")
   }, character(1L))

   out <- paste(names(x), "", x, sep = "\t")
   out <- paste(out, collapse = "\n")

   write(out, file = path, append = FALSE)
}

# Recommended to use R.utils::gzip on the resulting GMT file
