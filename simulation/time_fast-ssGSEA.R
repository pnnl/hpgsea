library(fast.ssgsea)

source("simulation/scripts/function-generate_data.R")

for_comparison <- FALSE
openblas <- grepl("openblas", sessionInfo()["BLAS"])

## Parameter combinations ----
if (for_comparison) {
   ## For comparison with ssGSEA-base
   param_list <- list("nGenes" = 1e4L,
                      "nSamples" = c(20L, 30L),
                      "minSetSize" = 10L,
                      "maxSetSize" = c(500L, 1000L),
                      "nSets" = c(1e3L, 2e3L),
                      "nperm" = c(1e3L, 2e3L),
                      "alpha" = c(0, 1))
} else {
   ## Extreme combinations to showcase speed of fast-ssGSEA
   param_list <- list("nGenes" = 1e4L,
                      "nSamples" = c(20L, 100L),
                      "minSetSize" = 10L,
                      "maxSetSize" = c(500L, 1000L),
                      "nSets" = c(1e3L, 1e4L, 1e5L),
                      "nperm" = c(1e3L, 1e4L),
                      "alpha" = c(0, 1))
}

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

      li <- with(row_i,
                 generate_data(nGenes, nSamples,
                               minSetSize, maxSetSize, nSets,
                               nperm, alpha))

      X_i <- li[["X"]]
      gene_sets_i <- li[["gene_sets"]]

      elapsed_time <- system.time({
         fast_ssgsea(X = X_i,
                     gene_sets = gene_sets_i,
                     alpha = row_i[["alpha"]],
                     nperm = row_i[["nperm"]],
                     batch_size = 1000L, # default
                     min_size = row_i[["minSetSize"]],
                     seed = 0L)
      })["elapsed"]

      message("  ", hms::as_hms(elapsed_time))

      comb_df_j[i, "elapsed_time"] <- elapsed_time
   }

   return(comb_df_j)
})

time_df <- do.call(what = rbind, args = time_df)

# Save results
if (openblas) {
   if (for_comparison) {
      ## First param_list, OpenBLAS
      saveRDS(
         object = time_df,
         file = file.path(
            "simulation",
            "data",
            "fast-ssGSEA_timing_results_OpenBLAS_for_comparison.rds"
         ),
         compress = TRUE,
         version = 3L
      )
   } else {
      ## Second param_list, OpenBLAS
      saveRDS(
         object = time_df,
         file = file.path("simulation",
                          "data",
                          "fast-ssGSEA_timing_results_OpenBLAS.rds"),
         compress = TRUE,
         version = 3L
      )
   }
} else {
   if (for_comparison) {
      ## First param_list, reference BLAS
      saveRDS(
         object = time_df,
         file = file.path("simulation",
                          "data",
                          "fast-ssGSEA_timing_results_BLAS_for_comparison.rds"),
         compress = TRUE,
         version = 3L
      )
   } else {
      ## Second param_list, reference BLAS
      saveRDS(
         object = time_df,
         file = file.path("simulation",
                          "data",
                          "fast-ssGSEA_timing_results_BLAS.rds"),
         compress = TRUE,
         version = 3L
      )
   }
}
