library(hpgsea)
library(dplyr)
library(cli)

source("manuscript/functions/function-create_dir.R")

data_dir <- create_data_dir()

gmt_file <- file.path(
  data_dir,
  "ptm.sig.db.all.flanking.human.v2.0.0.gmt.gz"
)

if (!file.exists(gmt_file)) {
  gmt_source <- paste0(
    "https://github.com/broadinstitute/ssGSEA2.0/blob/master/db/",
    "ptmsigdb/v2.0.0/all/ptm.sig.db.all.flanking.human.v2.0.0.gmt"
  )

  gzip_code <- sprintf(
    "R.utils::gzip(\"%s/%s\")",
    data_dir,
    basename(gmt_source)
  )

  error_msg <- paste0(
    "{gmt_file} does not exist. ",
    "Download the GMT file from {.url {gmt_source}}. ",
    "Then, compress the file with {.code {gzip_code}}."
  )

  cli_abort(message = error_msg)
}

x <- read_gmt(gmt_file)

# Remove flanking sequences that are both up and down in the same set
sets <- data.frame(
  set = rep(names(x), lengths(x)),
  elements = unlist(x)
) %>%
  mutate(
    elements_no_dir = substr(elements, 1L, 17L)
  ) %>%
  mutate(
    .by = c(
      set,
      elements_no_dir
    ),
    n = n()
  ) %>%
  filter(
    n == 1L
  ) %>%
  distinct(
    set,
    elements
  ) %>%
  summarise(
    .by = set,
    elements = list(elements)
  ) %>%
  tibble::deframe()

saveRDS(
  object = sets,
  file = file.path(data_dir, "ptm.sig.db.all.flanking.human.v2.0.0.rds"),
  compress = "xz"
)
