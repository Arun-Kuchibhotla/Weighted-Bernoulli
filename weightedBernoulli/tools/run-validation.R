#!/usr/bin/env Rscript
# Run the full base-R regression suite against the installed package.
# Usage: Rscript tools/run-validation.R [output_directory]
# Without an output directory, the report is placed in a temporary directory.

local({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 1L && args[1L] %in% c("-h", "--help")) {
    cat("Usage: Rscript tools/run-validation.R [output_directory]\n",
        "Install weightedBernoulli first. The optional directory receives ",
        "unified_validation.csv and randomized_validation.csv.\n", sep = "")
  } else {
    if (length(args) > 1L) {
      stop("Expected at most one argument: the validation output directory.")
    }
    script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                       value = TRUE)
    if (length(script_arg) != 1L) {
      stop("Run this script with Rscript tools/run-validation.R [output_directory].")
    }
    script_file <- normalizePath(sub("^--file=", "", script_arg),
                                 mustWork = TRUE)
    test_files <- file.path(dirname(dirname(script_file)), "tests",
                            c("regression.R", "randomized.R"))
    if (!all(file.exists(test_files))) {
      stop("Cannot find both regression suites beside the tools directory.")
    }
    if (!requireNamespace("weightedBernoulli", quietly = TRUE)) {
      stop("Install weightedBernoulli before running its package validation suite.")
    }
    output_dir <- if (length(args)) args[1L] else
      file.path(tempdir(), "weightedBernoulli-validation")
    if (!dir.exists(output_dir) &&
        !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Cannot create the validation output directory: ", output_dir)
    }
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    old_dir <- Sys.getenv("WEIGHTEDBERNOULLI_VALIDATION_DIR", unset = NA_character_)
    started <- proc.time()[["elapsed"]]
    tryCatch({
      Sys.setenv(WEIGHTEDBERNOULLI_VALIDATION_DIR = output_dir)
      for (test_file in test_files)
        sys.source(test_file, envir = new.env(parent = globalenv()))
    }, finally = {
      if (is.na(old_dir)) {
        Sys.unsetenv("WEIGHTEDBERNOULLI_VALIDATION_DIR")
      } else {
        Sys.setenv(WEIGHTEDBERNOULLI_VALIDATION_DIR = old_dir)
      }
    })
    cat(sprintf("Elapsed time: %.2f seconds.\n", proc.time()[["elapsed"]] - started))
    cat("Validation report: ", file.path(output_dir, "unified_validation.csv"),
        "\n", sep = "")
    cat("Randomized report: ", file.path(output_dir, "randomized_validation.csv"),
        "\n", sep = "")
  }
})
