suppressPackageStartupMessages({
  library(dplyr)
  library(hpgsea)

  source("manuscript/functions/function-ssgsea2.R")
  source("manuscript/functions/function-create_dir.R")
})

data_dir <- create_data_dir()

# See manuscript/prepare_PTMsigDB.R for details
sets <- file.path(data_dir, "ptm.sig.db.all.flanking.human.v2.0.0.rds") %>%
  readRDS()

flanking <- unlist(sets) %>%
  substr(1L, 17L) %>%
  unique()

n <- length(flanking)

set.seed(0L)
stats <- structure(rnorm(n), names = flanking)

# Benchmark parameters
nperm <- c(1e3L, 1e4L, 1e5L, 1e6L)
min_size <- 5L
alpha <- 1


# HPGSEA ----
time_hpgsea <- lapply(seq_len(3L), function(j) { # 3 replicates
  # Randomize order of runs
  set.seed(j)
  comb_df_j <- data.frame(
    nperm = sample(nperm),
    replicate = j
  )

  for (i in seq_len(nrow(comb_df_j))) {
    message(i)

    row_i <- comb_df_j[i, , drop = FALSE]

    invisible(gc())

    set.seed(0L)

    tic <- Sys.time()

    res <- hpgsea(
      stats = stats,
      gene_sets = sets,
      alpha = alpha,
      nperm = row_i[["nperm"]],
      min_size = min_size
    )

    toc <- Sys.time()

    elapsed_time <- as.numeric(
      difftime(toc, tic, units = "secs")
    )

    rm(list = "res")

    message("  ", hms::as_hms(elapsed_time))

    comb_df_j[i, "elapsed_time"] <- elapsed_time
  }

  return(comb_df_j)
})

time_hpgsea <- do.call(what = rbind, args = time_hpgsea)


## PTM-SEA (ssGSEA2.0) ----
nperm <- c(1e3L, 1e4L)

time_ssgsea2 <- lapply(seq_len(3L), function(j) { # 3 replicates
  comb_df_j <- data.frame(
    nperm = nperm,
    replicate = j
  )

  for (i in seq_len(nrow(comb_df_j))) {
    message(i)

    row_i <- comb_df_j[i, , drop = FALSE]

    invisible(gc())

    elapsed_time <- ssgsea2(
      stats = stats,
      gene_sets = sets,
      alpha = 1,
      nperm = row_i[["nperm"]],
      min_size = min_size,
      seed = 0L,
      time_only = TRUE
    )[[1L]]

    message("  ", hms::as_hms(elapsed_time))

    comb_df_j[i, "elapsed_time"] <- elapsed_time
  }

  return(comb_df_j)
})

time_ssgsea2 <- do.call(what = rbind, args = time_ssgsea2)

res <- list(
  "hpgsea" = time_hpgsea,
  "ssGSEA2.0" = time_ssgsea2
) %>%
  bind_rows(.id = "method")

saveRDS(
  object = res,
  file = file.path(data_dir, "time_PTM-SEA.rds")
)
