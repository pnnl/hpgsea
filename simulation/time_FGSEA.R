library(fgsea)

source("simulation/scripts/functions/function-generate_data.R")

# Parameter combinations ----
param_list <- list(
  "maxSetSize" = c(500L, 1000L),
  "nSets" = c(2e4L, 4e4L, 6e4L),
  "nperm" = c(1e4L, 1e5L, 1e6L)
)

n_genes <- 1e4L
min_size <- 10L
alpha <- 1

comb_df <- expand.grid(param_list)

# nproc = 0 for multithreading; nproc = 1 for single threading
for (nproc in c(1L, 0L)) {
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

      tic <- Sys.time()

      invisible({
        capture.output({
          res <- fgseaSimple(
            stats = stats_i,
            pathways = gene_sets_i,
            gseaParam = alpha,
            nperm = row_i[["nperm"]],
            minSize = min_size,
            maxSize = row_i[["maxSetSize"]],
            nproc = nproc
          )
        })
      })

      toc <- Sys.time()

      elapsed_time <- as.numeric(
        difftime(toc, tic, units = "secs")
      )

      message("  ", hms::as_hms(elapsed_time))

      comb_df_j[i, "elapsed_time"] <- elapsed_time

      rm(list = "res")
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
      sprintf("fgsea_times_nproc_%d.rds", nproc)
    )
  )

  Sys.sleep(30) # CPU cool down
}
