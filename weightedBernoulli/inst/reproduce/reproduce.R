# Reproduce the tables and figures for the ONE public quantile function.
# Run in a writable output directory, or use tools/reproduce.R [output_directory].
# Only standard R packages are required. No simulation is used.
local({
  if (!requireNamespace("weightedBernoulli", quietly = TRUE)) {
    stop("Install weightedBernoulli before running the comparison script.")
  }
  qweighted_bernoulli <- weightedBernoulli::qweighted_bernoulli
  dir.create("results", showWarnings = FALSE)
  dir.create("figures", showWarnings = FALSE)
  options(stringsAsFactors = FALSE, digits = 10)
  writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

  # Independent finite-distribution references, for this script only.
  lattice_reference <- function(w, pr, divisor = 1) {
    stopifnot(all(w >= 0), all(w == floor(w)))
    pr <- rep_len(pr, length(w))
    mass <- 1
    for (i in seq_along(w)) {
      mass <- c(mass * (1 - pr[i]), rep(0, w[i])) +
        c(rep(0, w[i]), mass * pr[i])
    }
    x <- seq.int(0, sum(w)) / divisor
    tail <- c(rev(cumsum(rev(mass[-1L]))), 0)
    list(x = x, mass = mass, tail = tail)
  }
  count_reference <- function() {
    k <- 0:300
    x <- as.vector(outer(k * sqrt(2), k * sqrt(3), "+"))
    mass <- as.vector(outer(stats::dbinom(k, 300, 0.2),
                            stats::dbinom(k, 300, 0.7), "*"))
    o <- order(x)
    x <- x[o]
    mass <- mass[o]
    list(x = x, mass = mass,
         tail = c(rev(cumsum(rev(mass[-1L]))), 0))
  }
  reference_quantile <- function(ref, alpha) {
    vapply(alpha, function(a) ref$x[which(ref$tail <= a)[1L]], numeric(1))
  }
  reference_tail <- function(ref, q) {
    vapply(q, function(x) sum(ref$mass[ref$x > x]), numeric(1))
  }
  short_method <- function(method) {
    map <- c(constant = "Constant", binomial = "Binomial",
      binomial_blocks = "Count classes", wilcoxon = "Ranks",
      lattice_dp = "Lattice DP", meet_in_the_middle = "MITM",
      grid = "Grid", chernoff = "Chernoff", hoeffding = "Hoeffding",
      cantelli = "Cantelli", bernstein = "Bernstein", binomial_p2 = "P2",
      bentkus_gaussian = "BD", bentkus_gaussian_k = "K-BD",
      pinelis = "Pinelis", pinelis_k = "K-Pinelis",
      symmetric_chebyshev = "Sym. Chebyshev", tomaszewski = "Tomaszewski",
      symmetry = "Symmetry", support = "Support", endpoint = "Top atom",
      lower_endpoint = "Lower endpoint", upper_endpoint = "Upper endpoint")
    z <- unname(map[method])
    z[grepl("^grouped_", method)] <- "Count envelope"
    z[grepl("^conditional_moments_", method)] <- "Conditional moments"
    z[grepl("^conditional_mgf_", method)] <- "Conditional MGF"
    z[is.na(z)] <- method[is.na(z)]
    z
  }
  family <- function(method) {
    z <- rep("Concentration", length(method))
    z[method == "grid"] <- "Grid"
    z[grepl("^grouped_", method)] <- "Count envelope"
    z[grepl("^conditional_moments_", method)] <- "Conditional moments"
    z[grepl("^conditional_mgf_", method)] <- "Conditional MGF"
    z
  }
  fixed <- function(x, digits = 3) {
    ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f"))
  }
  write_table <- function(path, columns, rows, align) {
    writeLines(c(paste0("\\begin{tabular}{", align, "}"), "\\toprule",
      paste0(paste(columns, collapse = " & "), " \\\\"), "\\midrule",
      vapply(rows, function(row) paste0(paste(row, collapse = " & "), " \\\\"),
             character(1)), "\\bottomrule", "\\end{tabular}"), path)
  }

  alpha <- c(0.05, 0.001, 1e-6)
  settings <- list(
    sqrt25 = list(label = "Square roots, 25", w = sqrt(1:25), pr = 0.5),
    sqrt26 = list(label = "Square roots, 26", w = sqrt(1:26), pr = 0.5),
    equal100 = list(label = "Equal, fair", w = rep(1, 100), pr = 0.5,
                    integer = rep(1, 100)),
    ranks100 = list(label = "Ranks, fair", w = 1:100, pr = 0.5, integer = 1:100),
    dominant40 = list(label = "60 + 39 unit weights", w = c(60, rep(1, 39)),
                      pr = 0.5, integer = c(60, rep(1, 39))),
    heterogeneous80 = list(label = "Lattice, unequal probabilities",
      w = rep(1:8, 10), pr = rep(seq(0.1, 0.8, length.out = 8), 10),
      integer = rep(1:8, 10)),
    repeated600 = list(label = "Two repeated real weights",
      w = c(rep(sqrt(2), 300), rep(sqrt(3), 300)),
      pr = c(rep(0.2, 300), rep(0.7, 300)), counts = TRUE),
    sqrt100 = list(label = "Square roots, 100", w = sqrt(1:100), pr = 0.5),
    near100 = list(label = "Near equal, fair", w = (1881 + 2 * (0:99)) / 1980,
                   pr = 0.5, integer = 1881 + 2 * (0:99), divisor = 1980),
    wide100 = list(label = "Spread weights, fair", w = (99 + 2 * (0:99)) / 198,
                   pr = 0.5, integer = 99 + 2 * (0:99), divisor = 198),
    near100_rare = list(label = "Near equal, success 0.1",
      w = (1881 + 2 * (0:99)) / 1980, pr = 0.1,
      integer = 1881 + 2 * (0:99), divisor = 1980),
    harmonic100 = list(label = "Harmonic, fair", w = 1 / (1:100), pr = 0.5),
    near1000 = list(label = "Near equal, 1000", w = seq(0.99, 1.01, length.out = 1000),
                    pr = 0.5),
    unequal_sqrt80 = list(label = "Square roots, unequal probabilities",
      w = sqrt(1:80), pr = seq(0.1, 0.8, length.out = 80))
  )
  comparisons <- bounds <- references <- fits <- list()
  for (id in names(settings)) {
    s <- settings[[id]]
    pr <- rep_len(s$pr, length(s$w))
    mu <- sum(s$w * pr)
    sd <- sqrt(sum(s$w^2 * pr * (1 - pr)))
    ref <- NULL
    if (!is.null(s$integer)) {
      divisor <- if (is.null(s$divisor)) 1 else s$divisor
      ref <- lattice_reference(s$integer, pr, divisor)
    }
    if (isTRUE(s$counts)) ref <- count_reference()
    if (!is.null(ref)) references[[id]] <- ref
    seconds <- system.time(fit <- qweighted_bernoulli(
      alpha, s$w, pr, lower.tail = FALSE, details = TRUE
    ))[["elapsed"]]
    fits[[id]] <- fit
    qref <- achieved <- rep(NA_real_, length(alpha))
    if (!is.null(ref)) {
      qref <- reference_quantile(ref, alpha)
      achieved <- reference_tail(ref, fit$quantile)
      stopifnot(all(fit$quantile >= qref - 1e-10 * max(1, sum(s$w))),
                all(achieved <= alpha * (1 + 1e-10) + 1e-15))
      exact <- fit$summary$exact
      stopifnot(all(abs(fit$quantile[exact] - qref[exact]) <=
                      1e-10 * max(1, sum(s$w))))
    }
    comparisons[[id]] <- data.frame(setting = id, label = s$label,
      n = length(s$w), alpha = alpha, mean = mu, sd = sd,
      reference_quantile = qref, quantile = fit$quantile,
      exact = fit$summary$exact, method = fit$summary$method,
      error_allowance = fit$summary$grid_error, actual_tail = achieved,
      seconds_for_three_levels = seconds)
    if (!is.null(fit$bounds)) {
      z <- fit$bounds
      z$setting <- id
      z$alpha <- alpha[z$request]
      z$family <- family(z$method)
      z$actual_tail <- if (is.null(ref)) NA_real_ else reference_tail(ref, z$quantile)
      ok <- z$available & is.finite(z$quantile)
      if (!is.null(ref)) {
        stopifnot(all(z$actual_tail[ok] <= z$tail_bound[ok] * (1 + 1e-10) + 1e-15))
      }
      for (j in seq_along(alpha)) {
        rows <- z$request == j & ok
        if (any(rows)) stopifnot(fit$quantile[j] == min(z$quantile[rows]))
      }
      bounds[[id]] <- z
    }
    cat("Completed", id, "in", seconds, "seconds.\n")
  }
  comparisons <- do.call(rbind, comparisons)
  bounds <- do.call(rbind, bounds)
  rownames(comparisons) <- rownames(bounds) <- NULL
  write.csv(comparisons, "results/comparisons.csv", row.names = FALSE)
  write.csv(bounds, "results/all_bounds.csv", row.names = FALSE)

  # Every result in the main table is the actual unified function output.
  main <- comparisons[comparisons$alpha == 0.05, ]
  rows <- lapply(seq_len(nrow(main)), function(i) {
    z <- main[i, ]
    c(z$label, z$n, fixed(z$reference_quantile), fixed(z$quantile),
      if (z$exact) "Yes" else "No", short_method(z$method))
  })
  write_table("results/main_table.tex",
    c("Setting", "$B$", "Reference", "Returned", "Exact", "Selected engine"),
    rows, "lrrrll")

  # Family minima are diagnostics from the same calls, not separate methods.
  families <- c("Concentration", "Grid", "Count envelope",
                "Conditional moments", "Conditional MGF")
  family_rows <- list()
  for (id in unique(bounds$setting)) for (j in seq_along(alpha)) {
    z <- bounds[bounds$setting == id & bounds$request == j & bounds$available, ]
    vals <- vapply(families, function(f) {
      q <- z$quantile[z$family == f]
      if (length(q)) min(q) else NA_real_
    }, numeric(1))
    selected <- comparisons[comparisons$setting == id & comparisons$alpha == alpha[j], ]
    family_rows[[length(family_rows) + 1L]] <- data.frame(
      setting = id, label = settings[[id]]$label, alpha = alpha[j],
      concentration = vals[1], grid = vals[2], envelope = vals[3],
      moments = vals[4], mgf = vals[5], returned = selected$quantile,
      method = selected$method, sd = selected$sd)
  }
  family_results <- do.call(rbind, family_rows)
  rownames(family_results) <- NULL
  write.csv(family_results, "results/family_comparisons.csv", row.names = FALSE)
  display <- family_results[family_results$alpha == 0.05, ]
  rows <- lapply(seq_len(nrow(display)), function(i) {
    z <- display[i, ]
    c(z$label, fixed(z$concentration, 3), fixed(z$grid, 3),
      fixed(z$envelope, 3), fixed(z$moments, 3), fixed(z$mgf, 3), fixed(z$returned, 3))
  })
  write_table("results/family_table.tex",
    c("Setting", "Concentr.", "Grid", "Counts", "Cond. mom.", "Cond. MGF", "Returned"),
    rows, "lrrrrrr")

  # Tail ratios at moderate and small levels, for computable reference laws.
  z <- comparisons[is.finite(comparisons$actual_tail) &
    comparisons$setting %in% c("equal100", "dominant40", "near100", "wide100", "near100_rare") &
    comparisons$alpha %in% c(0.05, 1e-6), ]
  rows <- lapply(seq_len(nrow(z)), function(i) {
    r <- z[i, ]
    a <- if (r$alpha == 1e-6) "$10^{-6}$" else "0.05"
    c(r$label, a, fixed(r$reference_quantile), fixed(r$quantile),
      fixed(r$actual_tail / r$alpha), short_method(r$method))
  })
  write_table("results/tail_table.tex",
    c("Setting", "$\\alpha$", "Reference", "Returned", "Actual tail/$\\alpha$", "Engine"),
    rows, "llrrrl")

  # Individual concentration candidates on two unstructured settings.
  methods <- c("chernoff", "hoeffding", "cantelli", "bernstein", "binomial_p2",
    "bentkus_gaussian", "bentkus_gaussian_k", "pinelis", "pinelis_k",
    "symmetric_chebyshev", "support")
  rows <- lapply(methods, function(method) {
    vals <- vapply(c("sqrt100", "harmonic100"), function(id) {
      z <- bounds[bounds$setting == id & bounds$request == 1 &
                  bounds$method == method & bounds$available, ]
      if (nrow(z)) min(z$quantile) else NA_real_
    }, numeric(1))
    c(short_method(method), fixed(vals[1], 4), fixed(vals[2], 4))
  })
  write_table("results/concentration_table.tex",
    c("Individual candidate", "Square roots, 100", "Harmonic, 100"), rows, "lrr")

  # A resource study still uses the same single function and automatic minimum.
  grid_study <- list()
  for (limit in c(512, 8192, 2^18)) {
    fit <- if (limit == 2^18) fits$sqrt100 else qweighted_bernoulli(
      0.05, sqrt(1:100), lower.tail = FALSE, details = TRUE,
      control = list(max_states = limit)
    )
    z <- subset(fit$bounds, request == 1 & method == "grid")
    grid_study[[length(grid_study) + 1L]] <- data.frame(
      max_states = limit, grid_quantile = z$quantile[1],
      grid_error = z$grid_error[1], returned = fit$quantile[1],
      method = fit$summary$method[1])
  }
  grid_study <- do.call(rbind, grid_study)
  write.csv(grid_study, "results/grid_resolution.csv", row.names = FALSE)
  rows <- lapply(seq_len(nrow(grid_study)), function(i) {
    z <- grid_study[i, ]
    c(format(z$max_states, scientific = FALSE, trim = TRUE),
      fixed(z$grid_quantile, 4), fixed(z$grid_error, 5),
      fixed(z$returned, 4), short_method(z$method))
  })
  write_table("results/grid_table.tex",
    c("State cap", "Grid candidate", "Grid allowance", "Returned", "Selected engine"),
    rows, "rrrrl")

  # Standard vector plots. Lines connect only the three evaluated levels.
  colors <- c(concentration = "#8D521A", grid = "#17846A", envelope = "#BB871F",
              moments = "#3978B3", mgf = "#9A529E", returned = "#182434")
  ltys <- c(concentration = 2, grid = 3, envelope = 4, moments = 5, mgf = 6, returned = 1)
  labels <- c(concentration = "Concentration min", grid = "Grid", envelope = "Count envelope",
              moments = "Conditional moments", mgf = "Conditional MGF", returned = "Returned")
  ids <- c("sqrt100", "wide100", "harmonic100", "unequal_sqrt80")
  grDevices::pdf("figures/bound_curves.pdf", width = 7.2, height = 6.3,
                family = "Helvetica", pointsize = 10)
  par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.8, 0.7), oma = c(0.2, 0, 1.8, 0),
      mgp = c(2.2, 0.65, 0))
  for (id in ids) {
    z <- family_results[family_results$setting == id, ]
    mu <- comparisons$mean[comparisons$setting == id][1]
    values <- (as.matrix(z[, names(colors)]) - mu) / z$sd
    plot(NA, xlim = range(-log10(z$alpha)), ylim = range(values, na.rm = TRUE),
      xlab = expression(-log[10](alpha)),
      ylab = expression((q-E[G])/sqrt(Var(G))),
      main = settings[[id]]$label, cex.main = 0.85)
    grid(col = "#E5E9EE", lty = 1)
    for (method in names(colors)) {
      if (all(is.na(values[, method]))) next
      lines(-log10(z$alpha), values[, method], col = colors[method],
            lty = ltys[method], lwd = if (method == "returned") 2 else 1.4,
            type = "b", pch = if (method == "returned") 16 else 1, cex = 0.55)
    }
    present <- names(colors)[colSums(is.finite(values)) > 0]
    legend("topleft", legend = labels[present], col = colors[present],
           lty = ltys[present], lwd = 1.5, cex = 0.64, bty = "n")
  }
  mtext("One function: candidate families and returned minimum", outer = TRUE,
        side = 3, line = 0.5, font = 2, cex = 1.05)
  grDevices::dev.off()

  grDevices::pdf("figures/improvement.pdf", width = 7.2, height = 4.7,
                family = "Helvetica", pointsize = 10)
  par(mar = c(4.3, 12, 2.9, 0.8), mgp = c(2.3, 0.65, 0))
  ids <- unique(family_results$setting)
  improvement <- vapply(alpha, function(a) {
    z <- family_results[family_results$alpha == a, ]
    z <- z[match(ids, z$setting), ]
    (z$concentration - z$returned) / z$sd
  }, numeric(length(ids)))
  yy <- rev(seq_along(ids))
  xrange <- range(0, improvement)
  plot(NA, xlim = xrange + c(-0.02, 0.05) * diff(xrange),
       ylim = c(0.4, length(ids) + 0.8), yaxt = "n", ylab = "",
       xlab = expression((q[concentration]-q[returned])/sqrt(Var(G))),
       main = "Reduction from the best concentration candidate", cex.main = 1.02)
  abline(v = pretty(xrange), col = "#E5E9EE")
  axis(2, at = yy, labels = vapply(settings[ids], function(s) s$label, character(1)),
       las = 1, cex.axis = 0.77, tick = FALSE)
  palette <- c("#17846A", "#3978B3", "#9A529E")
  for (j in seq_along(alpha)) points(improvement[, j], yy + (2 - j) * 0.15,
                                    pch = 14 + j, col = palette[j], cex = 0.8)
  legend("topright", legend = c("alpha = 0.05", "alpha = 0.001", "alpha = 1e-6"),
         pch = 15:17, col = palette, bty = "n", cex = 0.8)
  grDevices::dev.off()

  verified <- sum(bounds$available & is.finite(bounds$actual_tail))
  writeLines(c(
    paste("Main setting-level comparisons:", nrow(comparisons)),
    paste("Independent reference settings:", length(references)),
    paste("Exact returned setting-level results:", sum(comparisons$exact)),
    paste("Conservative returned setting-level results:", sum(!comparisons$exact)),
    paste("Individual candidate tail checks against reference laws:", verified),
    "Every returned conservative result equaled the minimum available candidate.",
    "Every checked strict tail satisfied bound*(1+1e-10)+1e-15.",
    "All exact returned values agreed with available independent reference quantiles.",
    "Rational lattice references describe the mathematically specified weights; arithmetic is double precision.",
    "Timings use one call with three probability levels; they are runtime descriptions, not native R benchmarks.",
    R.version.string, paste("Platform:", R.version$platform)
  ), "results/validation.txt")
  cat("Completed all unified comparisons and figures.\n")

  if (!grepl("emscripten", R.version$platform) && nzchar(Sys.which("gs"))) {
    for (file in c("figures/bound_curves.pdf", "figures/improvement.pdf")) {
      temporary <- paste0(file, ".embedded.pdf")
      ok <- try(grDevices::embedFonts(file, outfile = temporary), silent = TRUE)
      if (!inherits(ok, "try-error") && file.exists(temporary)) {
        file.copy(temporary, file, overwrite = TRUE)
        unlink(temporary)
      }
    }
  }
})
