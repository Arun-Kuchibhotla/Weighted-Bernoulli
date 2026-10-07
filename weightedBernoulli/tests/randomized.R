# Independent base-R validation of randomized cutoff distributions.
# R CMD check runs this file against the installed package.
# Every probability check uses an enumerated law or an explicit binomial law.

local({
# Independent finite-law references for randomized quantile validation.
# These enumerate bit vectors and evaluate the rejection events directly.
# They deliberately do not call a quantile, convolution, or package helper.

enumerate_bernoulli_reference <- function(weights, prob = 0.5) {
  prob <- rep_len(prob, length(weights))
  offset <- sum(weights[prob == 1])
  keep <- weights > 0 & prob > 0 & prob < 1
  weights <- weights[keep]
  prob <- prob[keep]
  n <- length(weights)
  if (n > 20L) stop("Reference enumeration is restricted to 20 bits.")
  if (!n) return(list(values = offset, mass = 1))
  bits <- vapply(seq_len(n), function(j) {
    bitwAnd(seq.int(0L, 2^n - 1L), bitwShiftL(1L, j - 1L)) != 0L
  }, logical(2^n))
  if (n == 1L) bits <- matrix(bits, ncol = 1L)
  values <- offset + as.vector(bits %*% weights)
  mass <- apply(bits, 1L, function(bit) prod(ifelse(bit, prob, 1 - prob)))
  keep <- mass > 0
  list(values = values[keep], mass = mass[keep])
}

reference_cutoff_distribution <- function(alpha, weights, prob = 0.5,
                                         two_sided = FALSE, law = NULL) {
  if (is.null(law)) law <- enumerate_bernoulli_reference(weights, prob)
  support <- sort(unique(law$values))
  if (two_sided) {
    total <- sum(weights)
    support <- c(-1, support[support <= total / 2])
    probability <- vapply(support, function(cutoff) {
      if (cutoff == -1) return(0)
      sum(law$mass[law$values <= cutoff | law$values >= total - cutoff])
    }, numeric(1))
    vapply(alpha, function(a) {
      k <- which(probability >= a)[1L]
      previous <- k - 1L
      c(lower = support[previous], upper = support[k],
        probability_lower = probability[previous],
        probability_upper = probability[k],
        tau = (a - probability[previous]) /
          (probability[k] - probability[previous]))
    }, numeric(5))
  } else {
    probability <- vapply(support, function(cutoff) {
      sum(law$mass[law$values > cutoff])
    }, numeric(1))
    vapply(alpha, function(a) {
      k <- which(probability <= a)[1L]
      previous <- k - 1L
      probability_lower <- if (previous) probability[previous] else 1
      c(lower = if (previous) support[previous] else -Inf,
        upper = support[k], probability_lower = probability_lower,
        probability_upper = probability[k],
        tau = (probability_lower - a) /
          (probability_lower - probability[k]))
    }, numeric(5))
  }
}

reference_randomized_rejection_probability <- function(
    lower, upper, tau, weights, prob = 0.5, two_sided = FALSE) {
  law <- enumerate_bernoulli_reference(weights, prob)
  reject <- if (two_sided) {
    total <- sum(weights)
    function(cutoff) {
      if (cutoff == -1) return(0)
      sum(law$mass[law$values <= cutoff | law$values >= total - cutoff])
    }
  } else {
    function(cutoff) sum(law$mass[law$values > cutoff])
  }
  tau * reject(upper) + (1 - tau) * reject(lower)
}

# This file is appended to references.R inside local({...}) for package tests.

qwb <- weightedBernoulli::qweighted_bernoulli
checks <- integer()
record <- function(category, n = 1L) {
  if (!category %in% names(checks)) checks[category] <<- 0L
  checks[category] <<- checks[category] + as.integer(n)
}
close <- function(x, y, tolerance = 2e-12) {
  isTRUE(length(x) == length(y) && all(x == y |
    (is.finite(x) & is.finite(y) &
     abs(x - y) <= tolerance * pmax(1, abs(x), abs(y)))))
}

assert_distribution <- function(fit, alpha, weights, prob = 0.5,
                                two_sided = FALSE, category,
                                reference = NULL) {
  if (is.null(reference)) reference <- reference_cutoff_distribution(
    alpha, weights, prob, two_sided)
  stopifnot(is.list(fit), is.data.frame(fit$summary),
            is.null(fit$kappa), is.null(fit$quantile), all(fit$exact),
            length(fit$c_plus) == length(alpha),
            close(fit$c_minus, reference["lower", ]),
            close(fit$c_plus, reference["upper", ]),
            close(fit$P_minus, reference["probability_lower", ]),
            close(fit$P_plus, reference["probability_upper", ]),
            close(fit$tau, reference["tau", ]),
            all(fit$tau >= 0 & fit$tau <= 1),
            close(fit$log_P_minus, log(fit$P_minus)),
            close(fit$log_P_plus, log(fit$P_plus)),
            close(fit$log_tau, log(fit$tau)),
            close(fit$log_one_minus_tau, log1p(-fit$tau)),
            close(fit$tau * fit$P_plus + (1 - fit$tau) * fit$P_minus,
                  alpha))
  if (length(weights) <= 20L) {
    actual <- vapply(seq_along(alpha), function(i) {
      reference_randomized_rejection_probability(
        fit$c_minus[i], fit$c_plus[i], fit$tau[i], weights, prob, two_sided)
    }, numeric(1))
    stopifnot(close(actual, alpha))
  }
  record(category, length(alpha))
  invisible(fit)
}

stopifnot(identical(getNamespaceExports("weightedBernoulli"),
                    "qweighted_bernoulli"),
          identical(formals(qwb)$randomized, FALSE),
          identical(formals(qwb)$two.sided, FALSE),
          identical(formals(qwb)$draw, TRUE))
record("public_randomized_API")

# Fair two-sided laws: every interior atom and between-atom level. Direct union
# evaluation includes the center once and treats the sentinel as an empty event.
for (n in 1:10) {
  w <- rep(1, n)
  d <- enumerate_bernoulli_reference(w)
  s <- sort(unique(d$values[d$values <= n / 2]))
  probability <- c(0, vapply(s, function(cutoff) {
    sum(d$mass[d$values <= cutoff | d$values >= n - cutoff])
  }, numeric(1)))
  alpha <- sort(unique(c(probability[probability > 0 & probability < 1],
                        (head(probability, -1L) + tail(probability, -1L)) / 2)))
  fit <- qwb(alpha, w, lower.tail = FALSE, randomized = TRUE,
             two.sided = TRUE, draw = FALSE)
  assert_distribution(fit, alpha, w, two_sided = TRUE,
                      category = "two_sided_all_jumps_and_midpoints")
}

for (w in list(c(1, 2), c(0.125, 0.25, 0.75, 1.5),
               c(0, 0.25, 0.25, 0.75, 1.5, 2), rep(0, 4))) {
  alpha <- c(0.03125, 0.125, 0.375, 0.625, 0.875)
  fit <- qwb(alpha, w, lower.tail = FALSE, randomized = TRUE,
             two.sided = TRUE, draw = FALSE)
  assert_distribution(fit, alpha, w, two_sided = TRUE,
                      category = "two_sided_gaps_fractional_and_zero")
}

# The original definition at its first jump and the central atom.
alpha <- c(0.05, 0.125, 0.625, 0.8)
fit <- qwb(alpha, rep(1, 4), lower.tail = FALSE, randomized = TRUE,
           two.sided = TRUE, draw = FALSE)
stopifnot(close(fit$c_minus, c(-1, -1, 0, 1)),
          close(fit$c_plus, c(0, 0, 1, 2)),
          close(fit$tau, c(0.4, 1, 1, 7 / 15)),
          fit$P_plus[4] == 1)
record("original_cutoff_example", length(alpha))

# One-sided unequal-probability laws: direct survival probabilities, including
# the predecessor below the minimum support point, at every atom and midpoint.
for (case in list(
    list(w = 1, pr = 0.125),
    list(w = c(1, 2, 4), pr = c(0.125, 0.5, 0.875)),
    list(w = c(0.25, 0.5, 1, 2, 4), pr = c(0.125, 0.25, 0.5, 0.75, 0.875)),
    list(w = c(0.25, 0.25, 1, 1.5), pr = c(0.25, 0.5, 0.75, 0.125)),
    list(w = c(0, 0.25, 0.5, 1, 2), pr = c(0.3, 0, 1, 0.25, 0.75)),
    list(w = c(0, 1, 2), pr = c(0.3, 1, 0)),
    list(w = c(0, 0), pr = c(0.125, 0.75)))) {
  d <- enumerate_bernoulli_reference(case$w, case$pr)
  probability <- sort(unique(c(0, 1, vapply(sort(unique(d$values)), function(q) {
    sum(d$mass[d$values > q])
  }, numeric(1)))))
  alpha <- sort(unique(c(probability[probability > 0 & probability < 1],
                        (head(probability, -1L) + tail(probability, -1L)) / 2)))
  fit <- qwb(alpha, case$w, case$pr, lower.tail = FALSE,
             randomized = TRUE, draw = FALSE)
  assert_distribution(fit, alpha, case$w, case$pr,
                      category = "one_sided_all_jumps_and_midpoints")
  # Complementing these dyadic levels is exact, so both argument conventions
  # must return the same distribution without a jump-boundary ambiguity.
  fit_lower <- qwb(1 - alpha, case$w, case$pr,
                   randomized = TRUE, draw = FALSE)
  stopifnot(close(fit_lower$c_plus, fit$c_plus),
            close(fit_lower$c_minus, fit$c_minus),
            close(fit_lower$tau, fit$tau))
  record("one_sided_input_conventions", length(alpha))
}

# Log input and the ordinary lower-tail convention in the two-sided mode.
w <- c(0.125, 0.25, 0.75, 1.5)
alpha <- c(0.03125, 0.375, 0.875)
fit <- qwb(1 - alpha, w, randomized = TRUE, two.sided = TRUE, draw = FALSE)
assert_distribution(fit, alpha, w, two_sided = TRUE,
                    category = "two_sided_lower_tail_input")
for (two_sided in c(FALSE, TRUE)) {
  alpha <- c(0.037, 0.213, 0.793)
  fit <- qwb(log(alpha), w, lower.tail = FALSE, log.p = TRUE,
             randomized = TRUE, two.sided = two_sided, draw = FALSE)
  assert_distribution(fit, alpha, w, two_sided = two_sided,
                      category = "log_probability_input")
}

# Long structured laws remain exact. The reference explicitly enumerates the
# independent binomial class counts rather than the Bernoulli bit vectors.
for (n in c(100, 101)) {
  w <- rep(0.25, n)
  d <- list(values = 0.25 * (0:n), mass = stats::dbinom(0:n, n, 0.5))
  alpha <- c(0.05, 0.2, 0.95, 0.999)
  reference <- reference_cutoff_distribution(alpha, w, two_sided = TRUE, law = d)
  fit <- qwb(alpha, w, lower.tail = FALSE, randomized = TRUE,
             two.sided = TRUE, draw = FALSE)
  assert_distribution(fit, alpha, w, two_sided = TRUE,
                      category = "long_fair_binomial", reference = reference)
}
w <- c(rep(0.25, 18), rep(0.75, 17))
pr <- c(rep(0.25, 18), rep(0.75, 17))
counts <- expand.grid(first = 0:18, second = 0:17)
d <- list(values = 0.25 * counts$first + 0.75 * counts$second,
          mass = stats::dbinom(counts$first, 18, 0.25) *
            stats::dbinom(counts$second, 17, 0.75))
alpha <- c(0.05, 0.2, 0.8)
reference <- reference_cutoff_distribution(alpha, w, pr, law = d)
fit <- qwb(alpha, w, pr, lower.tail = FALSE, randomized = TRUE, draw = FALSE)
assert_distribution(fit, alpha, w, pr, category = "long_repeated_classes",
                    reference = reference)

# A 25-bit fair binary-weight sum is uniform on integers 0,...,2^25-1.
# Its two-sided distribution is explicit; no full-state enumeration is needed.
w <- 2^(0:24)
alpha <- c(0.05, 0.5, 0.9)
c_plus <- ceiling(alpha * 2^24) - 1
P_plus <- (c_plus + 1) / 2^24
P_minus <- c_plus / 2^24
reference <- rbind(lower = c_plus - 1, upper = c_plus,
                   probability_lower = P_minus, probability_upper = P_plus,
                   tau = (alpha - P_minus) / (P_plus - P_minus))
fit <- qwb(alpha, w, lower.tail = FALSE, randomized = TRUE,
           two.sided = TRUE, draw = FALSE)
assert_distribution(fit, alpha, w, two_sided = TRUE,
                    category = "length25_randomized_exact", reference = reference)

# Log-scale calibration remains meaningful even when a displayed probability
# underflows to zero or the displayed tau rounds to one. Each law here has an
# explicit relevant endpoint atom, so the logarithmic mixture is independent.
for (case in list(
    list(level = -1000, w = 1, pr = 1e-100, two = FALSE, lower = 0, upper = 1),
    list(level = -1e-100, w = 1, pr = 0.5, two = FALSE, lower = -Inf, upper = 0),
    list(level = -1500, w = 2^(0:24), pr = 0.5, two = TRUE,
         lower = -1, upper = 0),
    list(level = -6000, w = 2^(0:24), pr = 1e-100, two = FALSE,
         lower = 2^25 - 2, upper = 2^25 - 1))) {
  fit <- qwb(case$level, case$w, case$pr, lower.tail = FALSE, log.p = TRUE,
             randomized = TRUE, two.sided = case$two, draw = FALSE)
  terms <- c(fit$log_tau + fit$log_P_plus,
             fit$log_one_minus_tau + fit$log_P_minus)
  log_actual <- max(terms) + log1p(exp(min(terms) - max(terms)))
  stopifnot(all(fit$exact), !fit$summary$degenerate,
            close(fit$c_minus, case$lower), close(fit$c_plus, case$upper),
            is.finite(log_actual), abs(log_actual / case$level - 1) < 1e-10)
  record("log_mixture_underflow_and_complement")
}

# Ordinary/default calls and distribution-only randomized queries consume no
# random numbers. A deterministic mixture consumes no random numbers either.
set.seed(1729)
before <- .Random.seed
invisible(qwb(c(0.5, 0.95), c(1, 2, 4)))
stopifnot(identical(.Random.seed, before))
invisible(qwb(c(0.05, 0.2), c(1, 2, 4), lower.tail = FALSE,
              randomized = TRUE, draw = FALSE))
stopifnot(identical(.Random.seed, before))
fit <- qwb(c(0.125, 0.625), rep(1, 4), lower.tail = FALSE,
           randomized = TRUE, two.sided = TRUE, draw = TRUE)
stopifnot(identical(.Random.seed, before), all(fit$tau == 1),
          identical(fit$kappa, fit$quantile), close(fit$kappa, fit$c_plus))
record("queries_and_degenerate_mixtures_preserve_RNG", 3L)

for (two_sided in c(FALSE, TRUE)) {
  alpha <- c(0.05, 0.125, 0.625, 0.8)
  set.seed(2718)
  fit <- qwb(alpha, rep(1, 4), lower.tail = FALSE, randomized = TRUE,
             two.sided = two_sided, draw = TRUE)
  after <- .Random.seed
  set.seed(2718)
  mixed <- which(fit$tau > 0 & fit$tau < 1)
  u <- stats::runif(length(mixed))
  expected <- fit$c_plus
  expected[mixed] <- ifelse(u <= fit$tau[mixed], fit$c_plus[mixed],
                           fit$c_minus[mixed])
  stopifnot(identical(.Random.seed, after), close(fit$kappa, expected),
            identical(fit$kappa, fit$quantile))
  set.seed(2718)
  again <- qwb(alpha, rep(1, 4), lower.tail = FALSE, randomized = TRUE,
               two.sided = two_sided, draw = TRUE)
  stopifnot(identical(again$kappa, fit$kappa))
  record("reproducible_independent_draws", length(alpha))
}

# Force the conservative branch on small laws that can be independently
# enumerated. Actual coverage and the direction of the stated certificate are
# checked directly; no approximate tail value is treated as an exact atom.
limited <- list(max_states = 1, max_dp_work = 1, max_mitm_work = 1,
                max_group_states = 0, max_count_states = 0, max_mgf_work = 0,
                max_p2_B = 0, max_p2_work = 0)
w <- c(1, 2, 5, 9, 17, 33, 65, 129)
for (two_sided in c(FALSE, TRUE)) {
  pr <- if (two_sided) 0.5 else (1:8) / 16
  alpha <- c(0.05, 0.2, 0.4)
  set.seed(57721)
  before <- .Random.seed
  fit <- qwb(alpha, w, pr, lower.tail = FALSE, max_length = 0,
             control = limited, randomized = TRUE, two.sided = two_sided,
             draw = TRUE)
  stopifnot(all(!fit$exact), close(fit$c_minus, fit$c_plus),
            all(fit$tau == 1), close(fit$kappa, fit$c_plus),
            identical(.Random.seed, before),
            all(fit$P_plus <= alpha + 2e-12),
            close(fit$P_plus, fit$P_minus))
  d <- enumerate_bernoulli_reference(w, pr)
  if (two_sided) stopifnot(all(fit$c_plus %in% c(-1, d$values)))
  actual <- vapply(seq_along(alpha), function(i) {
    reference_randomized_rejection_probability(
      fit$c_minus[i], fit$c_plus[i], fit$tau[i], w, pr, two_sided)
  }, numeric(1))
  stopifnot(all(actual <= alpha + 2e-12),
            all(actual <= fit$P_plus + 2e-12))
  record("conservative_actual_tail_and_RNG", length(alpha))
}

# Permit grouped numerical candidates while disabling exact lattice work.
# This case selects a grouped bound, so it verifies that every conservative
# candidate is evaluated at alpha/2 before forming a two-sided cutoff.
w <- 9:20
alpha <- c(0.05, 0.2)
fit <- qwb(alpha, w, lower.tail = FALSE, max_length = 0,
           control = list(max_states = 1, max_dp_work = 1, max_mitm_work = 1,
                          max_p2_B = 0, max_p2_work = 0),
           randomized = TRUE, two.sided = TRUE, draw = FALSE, details = TRUE)
selected <- fit$bounds[fit$bounds$selected, ]
grouped <- fit$bounds[fit$bounds$available &
                     grepl("grouped", fit$bounds$method), ]
stopifnot(all(!fit$exact), any(grepl("grouped", selected$method)),
          nrow(grouped) > 0,
          all(grouped$tail_bound <= alpha[grouped$request] / 2 + 2e-12))
actual <- vapply(seq_along(alpha), function(i) {
  reference_randomized_rejection_probability(
    fit$c_minus[i], fit$c_plus[i], fit$tau[i], w, two_sided = TRUE)
}, numeric(1))
stopifnot(all(actual <= alpha + 2e-12),
          all(actual <= fit$P_plus + 2e-12))
record("grouped_two_sided_alpha_half", length(alpha))

# Fairness and the interior-level requirements must be enforced explicitly.
throws <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
stopifnot(throws(qwb(0.05, c(1, 2), c(0.5, 0.4), lower.tail = FALSE,
                    randomized = TRUE, two.sided = TRUE)),
          throws(qwb(c(0, 1), rep(1, 3), lower.tail = FALSE,
                     randomized = TRUE, two.sided = TRUE)))
record("randomized_input_validation", 2L)
stopifnot(throws(qwb(0.5, c(2^54, 1, 2), c(1, 0.5, 0.5),
                    randomized = TRUE, draw = FALSE)))
record("collapsed_adjacent_support_refused")

report <- data.frame(category = names(checks), checks = unname(checks))
cat(sprintf("Randomized quantile validation: %d checks across %d categories passed.\n",
            sum(checks), length(checks)))
destination <- Sys.getenv("WEIGHTEDBERNOULLI_VALIDATION_DIR", unset = "")
if (nzchar(destination)) {
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(report, file.path(destination, "randomized_validation.csv"),
                   row.names = FALSE)
}

})
