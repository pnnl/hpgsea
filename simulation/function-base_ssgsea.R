# if (!require("package:cmapR", quietly = TRUE)) {
#    if (!require("package:BiocManager", quietly = TRUE))
#       install.packages("BiocManager", dependencies = TRUE)
#
#    BiocManager::install("cmapR", dependencies = TRUE)
# }
#
# if (!require("package:ssGSEA2", quietly = TRUE))
#    devtools::install_github("nicolerg/ssGSEA2", dependencies = TRUE)

library(ssGSEA2)
library(cmapR) # GCT, write_gct
library(dplyr)

source("simulation/scripts/function-writeGMT.R")

base_ssgsea <- function(X,
                        gene_sets,
                        alpha = 1,
                        nperm = 1000L,
                        min_size = 2L,
                        seed = NULL,
                        time_only = FALSE)
{
   set.seed(seed)

   # Create temporary directory to store input and results files
   temp_dir <- tempdir(check = TRUE)

   # Remove temporary directory on function exit
   on.exit(
      unlink(temp_dir, recursive = TRUE, force = TRUE)
   )

   gct_file <- file.path(temp_dir, "gct_file.gct")
   gmt_file <- file.path(temp_dir, "gmt_file.gmt")

   # Save list of gene sets to a GMT file
   writeGMT(gene_sets, gmt_file)

   # Save sample matrix to GCT file
   invisible(
      capture.output({ # do not print messages to console
         cmapR::write_gct(ds = cmapR::GCT(X),
                          ofile = gct_file,
                          precision = 8L, # affects ES calculation
                          appenddim = FALSE)
      })
   )

   args_list <- list(
      output.directory = temp_dir,
      input.ds = gct_file,
      gene.set.databases = gmt_file,
      output.prefix = "output",
      sample.norm.type = "none",
      weight = alpha,
      statistic = "area.under.RES",
      output.score.type = "NES",
      nperm = nperm,
      min.overlap = min_size,
      extended.output = FALSE,
      export.signat.gct = FALSE,
      param.file = FALSE,
      log.file = file.path(temp_dir, "log.log"),
      par = TRUE, # parallel processing
      spare.cores = 0L # use all threads
   )

   suppressMessages({
      invisible(
         capture.output({ # do not print messages to console
            elapsed_time <- system.time({
               ls <- do.call(what = ssGSEA2::run_ssGSEA2,
                             args = args_list)
            })["elapsed"]
         })
      )
   })

   out <- list("elapsed_time" = elapsed_time)

   if (!time_only) {
      # P-values
      pval_df <- file.path(temp_dir, "output-pvalues.gct") %>%
         read.delim(skip = 2L) %>%
         select(-starts_with("Signature")) %>%
         dplyr::rename(set = id) %>%
         tidyr::pivot_longer(cols = -set,
                             names_to = "sample",
                             values_to = "p_value")

      # Normalized enrichment scores
      NES_df <- file.path(temp_dir, "output-scores.gct") %>%
         read.delim(skip = 2L) %>%
         select(-starts_with("Signature")) %>%
         dplyr::rename(set = id) %>%
         tidyr::pivot_longer(cols = -set,
                             names_to = "sample",
                             values_to = "NES")

      # Enrichment scores
      ES_df <- ls %>%
         lapply(function(samples_i) {
            samples_i <- sapply(samples_i, function(xi) {
               xi[["ES"]]
            })

            out_i <- data.frame(ES = samples_i,
                                sample = colnames(X))
         }) %>%
         bind_rows(.id = "set") %>%
         select(set, sample, ES)

      df <- left_join(x = ES_df, y = NES_df,
                      by = c("sample", "set")) %>%
         left_join(y = pval_df, by = c("sample", "set"))

      out[["results"]] <- df
   }

   return(out)
}
