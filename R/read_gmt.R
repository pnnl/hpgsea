#' @title Construct a named list of gene sets from a GMT file
#'
#' @description Given a path to a Gene Matrix Transposed (GMT) file or GMT-like
#'   file, construct a named list of sets that can be used with
#'   \code{\link{hpgsea}}.
#'
#' @param file character; path to a GMT or GMT-like file. The file may be
#'   compressed.
#'
#' @seealso
#'   \href{https://docs.gsea-msigdb.org/#GSEA/Data_Formats/#gmt-gene-matrix-transposed-file-format-gmt}{GMT
#'   file format}
#'
#' @author Tyler Sagendorf
#'
#' @export read_gmt
#'
#' @examples
#' # HALLMARK gene sets from MSigDB
#' file <- system.file(
#'   "extdata",
#'   "h.all.v2025.1.Hs.symbols.gmt.gz",
#'   package = "hpgsea"
#' )
#'
#' gene_sets <- read_gmt(file)
#'
#' head(names(gene_sets)) # First 6 gene set names
#'
#' gene_sets[1] # first gene set
read_gmt <- function(file) {
  gmt <- readLines(con = file)
  gmt <- strsplit(gmt, split = "\t", fixed = TRUE)

  # Extract list of genes from each set. The first element is the name, and the
  # second is extra information about the set. The remaining elements are
  # genes.
  gene_sets <- lapply(gmt, function(x) x[3:length(x)])

  names(gene_sets) <- vapply(
    X = gmt,
    FUN = function(x) x[1L],
    FUN.VALUE = character(1L)
  )

  return(gene_sets)
}
