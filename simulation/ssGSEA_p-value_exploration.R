library(fast.ssgsea) # fast_ssgsea
library(dplyr)

source("simulation/scripts/function-generate_data.R")

# wrapper for ssGSEA2::run_ssGSEA2
source("simulation/scripts/function-base_ssgsea.R")


param_list <- list("nGenes" = 10000L,
                   "nSamples" = 20L,
                   "minSetSize" = 10L,
                   "maxSetSize" = 1000L,
                   "nSets" = 1000L,
                   "nperm" = 1000L,
                   "alpha" = c(0, 1))

comb_df <- expand.grid(param_list)


## Results ----

# Takes around 1 hour
tic <- Sys.time()

res_list <- lapply(seq_len(nrow(comb_df)), function(i) {
   message("i = ", i)

   row_i <- comb_df[i, , drop = FALSE]

   li <- with(
      data = row_i,
      expr = generate_data(nGenes = nGenes,
                           nSamples = nSamples,
                           minSetSize = minSetSize,
                           maxSetSize = maxSetSize,
                           nSets = nSets,
                           nperm = nperm,
                           alpha = alpha)
   )

   X_i <- li[["X"]]
   gene_sets_i <- li[["gene_sets"]]

   invisible(gc())

   fast_res <- fast_ssgsea(X = X_i,
                           gene_sets = gene_sets_i,
                           alpha = row_i[["alpha"]],
                           nperm = row_i[["nperm"]],
                           batch_size = 1e3L,
                           min_size = row_i[["minSetSize"]],
                           seed = 0L,
                           sort = FALSE) %>%
      select(-c(n_same_sign, n_as_extreme, adj_p_value))

   base_res <- base_ssgsea(X = X_i,
                           gene_sets = gene_sets_i,
                           alpha = row_i[["alpha"]],
                           nperm = row_i[["nperm"]],
                           min_size = row_i[["minSetSize"]],
                           seed = 0L,
                           time_only = FALSE)[["results"]]

   out <- list("fast" = fast_res,
               "base" = base_res) %>%
      # Modifications to reduce object size
      lapply(function(xi) {
         xi %>%
            mutate(across(.cols = c(sample, set),
                          .fns = ~ factor(.x, levels = unique(.x))),
                   across(.cols = c(ES, NES),
                          .fns = ~ signif(.x, digits = 6L)),
                   p_value = signif(p_value, digits = 4L))
      })

   return(out)
})

tok <- Sys.time()

message("Time elapsed: ", hms::as_hms(tok - tic))

names(res_list) <- comb_df$alpha

res_list <- lapply(res_list, function(li) {
   left_join(x = li[["base"]],
             y = li[["fast"]],
             by = c("sample", "set"),
             suffix = c(".base", ".fast"))
})

# with(res_list[["1"]],
#      plot(-log10(p_value.fast), -log10(p_value.base)))
# abline(a = 0, b = 1, col = "red")

# Save results for plotting
saveRDS(object = res_list,
        file = "simulation/data/ssGSEA_p-value_exploration_data.rds",
        compress = TRUE,
        version = 3L)
