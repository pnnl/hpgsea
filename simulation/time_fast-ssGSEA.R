library(fast.ssgsea)

source("simulation/scripts/function-generate_data.R")

# Parameter combinations ----
param_list <- list(
  "nGenes" = 1e4L,
  "nSamples" = 1L,
  "minSetSize" = 10L,
  "maxSetSize" = c(500L, 1000L),
  "nSets" = c(1e3L, 1e4L, 5e4L),
  "nperm" = c(1e4L, 1e5L, 1e6L),
  "alpha" = 1
)

comb_df <- expand.grid(param_list)

time_df <- lapply(seq_len(3L), function(j) { # 3 replicates
  # Randomize order of runs
  set.seed(j)
  comb_df_j <- comb_df[sample(seq_len(nrow(comb_df))), ]
  rownames(comb_df_j) <- NULL

  comb_df_j$replicate <- j

  for (i in seq_len(nrow(comb_df))) {
    message(i)

    row_i <- comb_df_j[i, , drop = FALSE]

    li <- with(
      row_i,
      generate_data(
        nGenes, nSamples,
        minSetSize, maxSetSize, nSets,
        nperm, alpha
      )
    )

    X_i <- li[["X"]]
    gene_sets_i <- li[["gene_sets"]]

    tic <- Sys.time()

    res <- fast_ssgsea(
      X = X_i,
      gene_sets = gene_sets_i,
      alpha = row_i[["alpha"]],
      nperm = row_i[["nperm"]],
      batch_size = 1e3L,
      min_size = row_i[["minSetSize"]],
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

time_df <- do.call(what = rbind, args = time_df)

# Save results
saveRDS(
  object = time_df,
  file = file.path(
    "simulation",
    "data",
    "fast-ssGSEA_timing_results.rds"
  )
)
