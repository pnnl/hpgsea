# This script simulates 10,000 sets of 10 to 500 genes to compare the runtime of
# HPGSEA to ssGSEA2.0.

suppressPackageStartupMessages({
  library(dplyr)
  library(hpgsea)

  source("simulation/scripts/functions/function-generate_data.R")
  source("simulation/scripts/functions/function-ssgsea2.R")
})


# Benchmark parameters
nperm <- c(1e3L, 1e4L, 1e5L, 1e6L)
min_size <- 10L
alpha <- 1

# Gene-level statistics and gene sets
li <- generate_data(
  nGenes = 1e4L,
  minSetSize = min_size,
  maxSetSize = 500L,
  nSets = 1e4L
)

stats <- li[["stats"]]
gene_sets <- li[["gene_sets"]]

# HPGSEA ----
time_hpgsea <- lapply(seq_len(3L), function(j) { # 3 replicates
  set.seed(j)
  comb_df_j <- data.frame(
    nperm = sample(nperm),
    replicate = j
  )

  for (i in seq_len(nrow(comb_df_j))) {
    message(i)

    row_i <- comb_df_j[i, , drop = FALSE]

    invisible(gc())

    tic <- Sys.time()

    res <- hpgsea(
      stats = stats,
      gene_sets = gene_sets,
      alpha = alpha,
      nperm = row_i[["nperm"]],
      min_size = min_size,
      seed = 0L
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

# ssGSEA2.0 ----
nperm <- c(1e3L, 1e4L)

# Results are continuously saved to a temporary file in case the analysis gets
# interrupted (e.g., power outage).
temp_file <- "simulation/data/TEMP_ssGSEA2.0_times.rds"

if (file.exists(temp_file)) {
  file.remove(temp_file)
}

time_vec <- numeric(0L)

time_ssgsea2 <- lapply(seq_len(3L), function(j) {
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
      gene_sets = gene_sets,
      alpha = 1,
      nperm = row_i[["nperm"]],
      min_size = min_size,
      seed = 0L,
      time_only = TRUE
    )[[1L]]

    time_vec <<- c(time_vec, elapsed_time)

    saveRDS(time_vec, file = temp_file)

    message("  ", hms::as_hms(elapsed_time))

    comb_df_j[i, "elapsed_time"] <- elapsed_time
  }

  return(comb_df_j)
})

time_ssgsea2 <- do.call(what = rbind, args = time_ssgsea2)

res <- list(
  "HPGSEA" = time_hpgsea,
  "ssGSEA2.0" = time_ssgsea2
) %>%
  bind_rows(.id = "method")

saveRDS(
  object = res,
  file = "simulation/data/ssGSEA2.0_times.rds"
)

if (file.exists(temp_file)) {
  file.remove(temp_file)
}
