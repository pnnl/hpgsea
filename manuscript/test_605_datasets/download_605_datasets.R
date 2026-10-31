# Fails if connected to a VPN

suppressPackageStartupMessages({
  library(rvest)
  library(dplyr)
  library(purrr)
  library(glue)
  library(data.table)
  library(float)
  library(R.utils)

  source("manuscript/functions/function-create_dir.R")
})

data_dir <- create_data_dir()

out_dir <- file.path(
  data_dir,
  "test_605_datasets"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

gmt_files <- file.path(
  out_dir,
  c(
    "c5.go.v2026.1.Hs.symbols.gmt",
    "m5.go.v2026.1.Mm.symbols.gmt"
  )
)

for (gmt_file in gmt_files) {
  if (!file.exists(paste0(gmt_file, ".gz"))) {
    root <- sub(".*go\\.v(.*)\\.symbols.*", "\\1", gmt_file)

    dest_file <- paste(out_dir, basename(gmt_file), sep = "/")

    download.file(
      url = paste(
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release",
        root,
        basename(gmt_file),
        sep = "/"
      ),
      destfile = dest_file,
      method = "libcurl"
    )

    gzip(dest_file)

    gmt <- read_gmt(paste0(dest_file, ".gz"))

    saveRDS(
      object = gmt,
      file = sub(".gmt", ".rds", dest_file, fixed = TRUE),
      compress = "xz"
    )
  }
}

get_file_names <- function(html) {
  files <- html_text2(html)
  files <- strsplit(files, split = "\t|\n")[[1L]]
  files <- files[-seq_len(14L)] # remove header info

  files <- stringr::str_trim(files) # remove whitespace
  files <- setdiff(files, "")

  # Remove junk
  keep <- !grepl("^\\d+K$", files) &
    !grepl("^2020-\\d{2}-\\d{2} .*", files) &
    !grepl("Ubuntu", files)

  files <- files[keep]

  return(files)
}

root_dir <- "https://ctlab.itmo.ru/files/software/fgsea/geo_ranks"

# This takes about 30 minutes
for (organism in c("mouse", "human")) {

  out_file <- file.path(
    out_dir,
    glue("datasets_{organism}.rds")
  )

  if (!file.exists(out_file)) {
    organism_root_dir <- file.path(
      root_dir,
      ifelse(organism == "human", "hsa", "mmu")
    )

    html_data <- read_html(organism_root_dir)

    files <- get_file_names(html_data)

    if (organism == "mouse") {
      # Fix one file name
      files[grepl("%", files)] <-
        "GSE27975_GPL1261.protocol.normoxia_21%25_O2.hypoxia_1%25_O2.rnk"
    }

    files <- file.path(organism_root_dir, files)

    stats_list <- map(
      .x = files,
      .f = function(f) {
        df <- fread(
          file = f,
          sep = "\t",
          showProgress = FALSE,
          verbose = FALSE
        )

        out <- list(
          "Symbol" = df$Symbol,
          "t" = as.float(df$t) # saves about 11 MB of disk space per organism
        )

        return(out)
      },
      .progress = list(
        name = glue("Downloading {organism} datasets"),
        clear = FALSE,
        show_after = 0
      )
    )

    names(stats_list) <- basename(files)

    saveRDS(
      object = stats_list,
      file = out_file,
      compress = "xz"
    )
  }

}
