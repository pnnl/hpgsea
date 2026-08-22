# Fails if connected to a VPN

suppressPackageStartupMessages({
  library(rvest)
  library(dplyr)
  library(purrr)
  library(glue)
  library(data.table)
  library(float)
})

out_dir <- file.path(
  "simulation",
  "data",
  "test_605_datasets"
)

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
    file = file.path(
      out_dir,
      glue("datasets_{organism}.rds")
    ),
    compress = "xz"
  )

}
