# if (!require("ssGSEA2", quietly = TRUE))
#    devtools::install_github("nicolerg/ssGSEA2")

library(cmapR) # GCT

source("simulation/scripts/function-generate_data.R")
source("simulation/scripts/function-base_ssgsea.R")


# Parameter combinations for comparison with ssGSEA-fast
param_list <- list("nGenes" = 10000L,
                   "nSamples" = c(20L, 30L),
                   "minSetSize" = 10L,
                   "maxSetSize" = c(500L, 1000L),
                   "nSets" = c(1000L, 2000L),
                   "nperm" = c(1000L, 2000L),
                   "alpha" = c(0, 1))

comb_df <- expand.grid(param_list)

time_df <- lapply(seq_len(1L), function(j) { # 1 replicate due to long runtimes
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

      elapsed_time <- base_ssgsea(X = X_i,
                                  gene_sets = gene_sets_i,
                                  alpha = row_i[["alpha"]],
                                  nperm = row_i[["nperm"]],
                                  min_size = row_i[["minSetSize"]],
                                  seed = 0L,
                                  time_only = TRUE)[["elapsed_time"]]

      message("  ", hms::as_hms(elapsed_time))

      # Save intermediate results, just in case
      write(x = elapsed_time,
            file = "simulation/data/base-ssGSEA_elapsed_time.txt",
            append = TRUE,
            sep = "\n")

      comb_df_j[i, "elapsed_time"] <- elapsed_time
   }

   return(comb_df_j)
})

time_df <- do.call(what = rbind, args = time_df)

# Save full results
saveRDS(
   object = time_df,
   file = file.path(
      "simulation",
      "data",
      "base-ssGSEA_timing_results_for_comparison.rds"
   ),
   compress = TRUE,
   version = 3L
)

# No longer needed
file.remove("simulation/data/base-ssGSEA_elapsed_time.txt")
