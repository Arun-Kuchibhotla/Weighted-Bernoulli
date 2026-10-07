# Base-R regression tests for the installed package.
# R CMD check runs this file automatically. To run it directly after installing
# the package, use Rscript tests/regression.R from the repository root.
# Set WEIGHTEDBERNOULLI_VALIDATION_DIR to write a CSV report; otherwise these
# tests leave the source tree and the installed package unchanged.

local({
qwb <- weightedBernoulli::qweighted_bernoulli
stopifnot(identical(getNamespaceExports("weightedBernoulli"),
                    "qweighted_bernoulli"),
          identical(as.numeric(formals(qwb)$max_length), 25),
          !"method" %in% names(formals(qwb)))

checks <- integer()
record <- function(category, number) {
  if (!category %in% names(checks)) checks[category] <<- 0L
  checks[category] <<- checks[category] + as.integer(number)
  invisible(NULL)
}
record("single_public_function_and_default", 2L)

log_add <- function(x, y) {
  z <- max(x, y)
  if (!is.finite(z)) return(z)
  z + log1p(exp(-abs(x - y)))
}
log_complement <- function(x) {
  ifelse(x < -log(2), log1p(-exp(x)), log(-expm1(x)))
}

# Independent exact reference for binary weights 1,2,4,...: compare binary
# digits from largest to smallest and sum the first-different-digit events.
# This reference never enumerates subsets and does not use a convolution.
binary_log_tail <- function(k, pr, upper = TRUE) {
  equal <- 0
  result <- -Inf
  for (j in rev(seq_along(pr))) {
    bit <- (floor(k / 2^(j - 1L)) %% 2) == 1
    if (upper && !bit) result <- log_add(result, equal + log(pr[j]))
    if (!upper && bit) result <- log_add(result, equal + log1p(-pr[j]))
    equal <- equal + if (bit) log(pr[j]) else log1p(-pr[j])
  }
  if (!upper) result <- log_add(result, equal)
  min(0, result)
}
binary_reference <- function(p, pr, lower.tail = TRUE, log.p = FALSE) {
  lp <- if (log.p) p else log(p)
  other <- if (log.p) log_complement(p) else log1p(-p)
  ll <- if (lower.tail) lp else other
  lu <- if (lower.tail) other else lp
  vapply(seq_along(p), function(i) {
    lo <- 0
    hi <- 2^length(pr) - 1
    use_upper <- ll[i] > -log(2)
    while (lo < hi) {
      mid <- floor(lo + (hi - lo) / 2)
      ok <- if (use_upper) binary_log_tail(mid, pr, TRUE) <= lu[i] else
        binary_log_tail(mid, pr, FALSE) >= ll[i]
      if (ok) hi <- mid else lo <- mid + 1
    }
    lo
  }, numeric(1))
}

assert_exact <- function(fit, reference, category, tolerance = 0) {
  stopifnot(all(fit$summary$exact), length(fit$quantile) == length(reference),
            all(is.finite(fit$quantile)),
            all(abs(fit$quantile - reference) <= tolerance))
  record(category, length(reference))
}

# Every interior fair-binomial CDF atom up to the promised length cutoff.
for (n in 2:25) {
  mass <- 1
  for (i in seq_len(n)) mass <- c(mass, 0) + c(0, mass)
  p <- head(cumsum(mass * 2^-n), -1L)
  fit <- qwb(p, rep(1, n), details = TRUE)
  assert_exact(fit, seq.int(0, n - 1L), "fair_binomial_CDF_atoms")
}

# Length 25: the dense lattice is too large, so this exercises the short-vector
# exact guarantee through two half arrays (4096 and 8192 entries), not 2^25
# full states. It also checks both sides of the median CDF atom.
w <- 2^(0:24)
p <- c(0.01, 0.10, 0.5 - 2^-53, 0.5, 0.5 + 2^-53, 0.95, 0.999)
ref <- ceiling(p * 2^25) - 1
fit <- qwb(p, w, details = TRUE)
assert_exact(fit, ref, "length25_fair_exact")
stopifnot(all(fit$summary$method == "meet_in_the_middle"))
record("length25_route", 1L)

pr <- seq(0.03, 0.97, length.out = 25)
p <- c(0.037, 0.213, 0.5, 0.793, 0.95, 0.997)
fit <- qwb(p, w, prob = pr, details = TRUE)
assert_exact(fit, binary_reference(p, pr), "length25_unequal_exact")

# Restrictive conservative work limits cannot change the default short-vector
# promise to an approximate answer.
fit <- qwb(c(0.5, 0.95), w, details = TRUE,
  control = list(max_states = 1, max_dp_work = 1, max_mitm_work = 1,
                 max_iterations = 1, max_group_states = 0))
assert_exact(fit, ceiling(c(0.5, 0.95) * 2^25) - 1,
             "length25_exact_despite_small_controls")

# Very small success probabilities expose both endpoint cancellation and
# underflow of products of Bernoulli probabilities.
fit <- qwb(c(1e-99, 1e-200), weights = 1, prob = 1e-100,
           lower.tail = FALSE, details = TRUE)
assert_exact(fit, c(0, 1), "tiny_probability_endpoint")
pr <- rep(1e-100, 25)
p <- c(1e-205, 1e-310)
fit <- qwb(p, w, prob = pr, lower.tail = FALSE, details = TRUE)
assert_exact(fit, binary_reference(p, pr, FALSE), "length25_tiny_upper_tail")
fit <- qwb(-1000, w, prob = pr, lower.tail = FALSE, log.p = TRUE,
           details = TRUE)
assert_exact(fit, binary_reference(-1000, pr, FALSE, TRUE),
             "length25_log_upper_tail")
pr <- rep(2^-1074, 25)
fit <- qwb(-1500, w, prob = pr, lower.tail = FALSE, log.p = TRUE,
           details = TRUE)
assert_exact(fit, binary_reference(-1500, pr, FALSE, TRUE),
             "length25_subnormal_probabilities")
pr <- rep(1 - 2^-40, 25)
fit <- qwb(-1000, w, prob = pr, log.p = TRUE, details = TRUE)
assert_exact(fit, binary_reference(-1000, pr, TRUE, TRUE),
             "length25_log_lower_tail")

# Common power-of-two rescaling must preserve the exact quantile.
for (scale in c(2^-500, 2^500)) {
  p <- c(0.25, 0.95)
  fit <- qwb(p, w * scale, details = TRUE)
  assert_exact(fit, (ceiling(p * 2^25) - 1) * scale,
               "exact_binary_rescaling")
}

# Independent exhaustive reference for small arbitrary weights; survival
# probabilities are accumulated directly, rather than subtracting a CDF.
enumerate <- function(w, pr) {
  values <- 0
  mass <- 1
  pr <- rep_len(pr, length(w))
  for (j in seq_along(w)) {
    values <- c(values, values + w[j])
    mass <- c(mass * (1 - pr[j]), mass * pr[j])
  }
  at <- order(values)
  values <- values[at]
  mass <- mass[at]
  list(values = values, mass = mass, cdf = cumsum(mass),
       survival = c(rev(cumsum(rev(mass[-1L]))), 0))
}
reference_quantile <- function(d, p) {
  vapply(p, function(x) {
    k <- if (x <= 0.5) which(d$cdf >= x)[1L] else
      which(d$survival <= 1 - x)[1L]
    d$values[k]
  }, numeric(1))
}
w <- sqrt(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29))
d <- enumerate(w, 0.5)
p0 <- seq.int(8, 1016, by = 16) / 1024
p <- sort(unique(c(p0 - 2^-53, p0, p0 + 2^-53)))
fit <- qwb(p, w, details = TRUE)
assert_exact(fit, reference_quantile(d, p), "arbitrary_weight_CDF_atoms",
             tolerance = 1e-12 * sum(w))

# Recognizable exact distributions above length 25.
p <- c(0.025, 0.5, 0.95)
assert_exact(qwb(p, rep(1, 1000), prob = 0.2, details = TRUE),
             stats::qbinom(p, 1000, 0.2), "large_binomial_structure")
assert_exact(qwb(p, 1:100, details = TRUE), stats::qsignrank(p, 100),
             "large_rank_structure")

# Three repeated integer weights/probabilities: enumerate binomial counts as
# a reference to the public lattice or block-count implementation.
count_reference <- function(a, size, success) {
  v <- 0
  m <- 1
  for (j in seq_along(a)) {
    k <- seq.int(0, size[j])
    old <- length(v)
    v <- as.vector(outer(v, a[j] * k, "+"))
    m <- rep(m, length(k)) * rep(stats::dbinom(k, size[j], success[j]), each = old)
  }
  at <- order(v)
  v <- v[at]
  m <- m[at]
  list(values = v, mass = m, cdf = cumsum(m),
       survival = c(rev(cumsum(rev(m[-1L]))), 0))
}
a <- c(2, 6, 10)
size <- c(20, 20, 20)
success <- c(0.2, 0.5, 0.8)
w <- rep(a, size)
pr <- rep(success, size)
d <- count_reference(a, size, success)
p <- c(0.05, 0.5, 0.95)
assert_exact(qwb(p, w, prob = pr, details = TRUE),
             reference_quantile(d, p), "large_lattice_structure")

# Irrational repeated coefficients with different success-probability classes.
a <- c(sqrt(2), pi)
size <- c(60, 40)
success <- c(0.2, 0.7)
w <- rep(a, size)
pr <- rep(success, size)
d <- count_reference(a, size, success)
p <- c(0.05, 0.5, 0.95)
fit <- qwb(p, w, prob = pr, details = TRUE)
assert_exact(fit, reference_quantile(d, p), "large_nonlattice_count_structure",
             tolerance = 1e-12 * sum(w))
stopifnot(all(fit$summary$method == "binomial_blocks"))
record("binomial_block_route", 1L)

# Deterministic/zero terms, endpoint conventions, names, and empty requests.
fit <- qwb(c(0, 0.5, 0.5 + 2^-53, 1), c(0, 2, 5, 9),
           prob = c(0.3, 0, 1, 0.5), details = TRUE)
assert_exact(fit, c(5, 5, 14, 14), "endpoints_and_deterministic_terms")
assert_exact(qwb(c(0, 0.5, 1), c(2, 5), prob = c(0, 1), details = TRUE),
             c(5, 5, 5), "fully_deterministic")
fit <- qwb(0.95, c(rep(0, 40), sqrt(c(2, 3, 5))), details = TRUE)
stopifnot(fit$summary$exact, fit$diagnostics$active_length == 3L)
record("active_length_reduction", 1L)
named_p <- c(low = 0.1, high = 0.9)
stopifnot(identical(names(qwb(named_p, 1:5)), names(named_p)),
          identical(qwb(numeric(), 1:5), numeric()))
record("names_and_empty_requests", 2L)

check_minimum <- function(fit, category) {
  stopifnot(!is.null(fit$bounds))
  for (i in seq_along(fit$quantile)) {
    b <- fit$bounds[fit$bounds$request == i & fit$bounds$available, ]
    stopifnot(nrow(b) > 0L, all(is.finite(b$quantile)),
              abs(fit$quantile[i] - min(b$quantile)) <=
                8 * .Machine$double.eps * max(1, abs(fit$quantile[i])),
              sum(b$selected) == 1L,
              fit$summary$method[i] == b$method[b$selected])
  }
  record(category, length(fit$quantile))
}

# Length 26, genuinely nonlattice weights: dispatch directly to the complete
# applicable bound set, and return its reported minimum.
w <- sqrt(seq_len(26))
fit <- qwb(c(0.95, 0.999), w, details = TRUE)
stopifnot(!any(fit$summary$exact), fit$diagnostics$active_length == 26L,
          is.null(fit$diagnostics$exact_attempts$meet_in_the_middle))
check_minimum(fit, "length26_minimum_of_candidates")

# Independently check every available strict-tail candidate against exhaustive
# small distributions. max_length=0 deliberately exercises the bound branch.
w <- sqrt(c(2, 3, 5, 7, 11, 13, 17, 19, 23))
control <- list(max_states = 4096, max_dp_work = 1000000,
                max_group_states = 128, max_count_states = 256,
                max_mgf_work = 1000000)
for (pr in list(0.5, 0.1, seq(0.1, 0.9, length.out = length(w)))) {
  d <- enumerate(w, pr)
  fit <- qwb(c(0.90, 0.99), w, prob = pr, max_length = 0,
             details = TRUE, control = control)
  check_minimum(fit, "forced_bound_minimum")
  b <- fit$bounds[fit$bounds$available, ]
  actual <- vapply(b$quantile, function(q) sum(d$mass[d$values > q]), numeric(1))
  stopifnot(all(actual <= b$tail_bound * (1 + 1e-10) + 1e-15))
  record("exhaustive_candidate_tail", length(actual))
}

validation <- data.frame(check = names(checks), count = unname(checks),
                         passed = TRUE, stringsAsFactors = FALSE)
output_dir <- Sys.getenv("WEIGHTEDBERNOULLI_VALIDATION_DIR", unset = "")
if (nzchar(output_dir)) {
  if (!dir.exists(output_dir) &&
      !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create the requested validation output directory: ", output_dir)
  }
  write.csv(validation, file.path(output_dir, "unified_validation.csv"),
            row.names = FALSE)
}
cat("Installed-package quantile validation passed:", sum(checks),
    "checks across", length(checks), "categories.\n")
print(validation, row.names = FALSE)
invisible(validation)
})
