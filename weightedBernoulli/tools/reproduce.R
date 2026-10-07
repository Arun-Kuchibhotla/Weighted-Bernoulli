#!/usr/bin/env Rscript
# Run every comparison and figure using the installed package.
# Usage: Rscript tools/reproduce.R [output_directory]
local({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Supply at most one output directory.")
  if (!requireNamespace("weightedBernoulli", quietly = TRUE)) {
    stop("Install weightedBernoulli before running the comparisons.")
  }
  output <- if (length(args)) args[1L] else "reproduced-results"
  if (!dir.exists(output) &&
      !dir.create(output, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create the output directory: ", output)
  }
  output <- normalizePath(output, mustWork = TRUE)
  script <- system.file("reproduce", "reproduce.R", package = "weightedBernoulli",
                        mustWork = TRUE)
  old <- setwd(output)
  tryCatch(sys.source(script, envir = new.env(parent = globalenv())),
           finally = setwd(old))
  cat("Comparison data and figures were written to ", output, ".\n", sep = "")
})
