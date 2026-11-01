library(hpgsea)

# To generate the "no_SIMD" timings, comment out the "#define USE_SIMD" line in
# src/Rcpp_functions.cpp. Then, build and install hpgsea and set USE_SIMD to
# FALSE.
USE_SIMD <- TRUE

source("manuscript/functions/function-generate_data.R")
source("manuscript/functions/function-create_dir.R")

data_dir <- create_data_dir()

# Parameter combinations ----
param_list <- list(
  "maxSetSize" = c(500L, 1000L),
  "nSets" = c(2e4L, 4e4L, 6e4L)
)

n_genes <- 1e4L
min_size <- 10L
alpha <- 1

comb_df <- expand.grid(param_list)

for (nthreads in c(1L, 0L)) {
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
          nGenes = n_genes,
          minSetSize = min_size,
          maxSetSize,
          nSets
        )
      )

      stats_i <- li[["stats"]]
      gene_sets_i <- li[["gene_sets"]]

      invisible(gc())

      set.seed(0L)

      tic <- Sys.time()

      res <- hpgsea(
        stats = stats_i,
        gene_sets = gene_sets_i,
        alpha = alpha,
        nperm = 1e6L,
        min_size = min_size,
        max_size = row_i[["maxSetSize"]],
        nthreads = nthreads
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

  if (!dir.exists(data_dir)) {
    invisible(dir.create(data_dir))
  }

  # Save results
  saveRDS(
    object = time_df,
    file = file.path(
      data_dir,
      sprintf("time_hpgsea_nthreads_%s%s.rds", nthreads,
              ifelse(USE_SIMD, "", "_no_SIMD"))
    )
  )
}
