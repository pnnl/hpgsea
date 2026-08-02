# The GMT file for PTMsigDB (version 2.0.0) was downloaded from
# https://github.com/broadinstitute/ssGSEA2.0/tree/master/db/ptmsigdb/v2.0.0/all

library(hpgsea)
library(dplyr)

# This file is 2.5 MB, so it is not tracked by git
x <- read_gmt("simulation/data/ptm.sig.db.all.flanking.human.v2.0.0.gmt")

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
  file = "simulation/data/ptm.sig.db.all.flanking.human.v2.0.0.rds",
  compress = "xz"
)
