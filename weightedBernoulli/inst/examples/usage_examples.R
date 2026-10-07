# Executable examples for the single public quantile function.
# Load the installed package, then run this script in a writable directory.
library(weightedBernoulli)

# Original fair two-sided cutoff: compute the mixing distribution first.
set.seed(20261007)
seed_before <- .Random.seed
cut <- qweighted_bernoulli(
  0.05, weights = 1:20, lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE, draw = FALSE
)
stopifnot(
  identical(seed_before, .Random.seed),
  is.null(cut$kappa), is.null(cut$quantile), cut$exact,
  cut$c_minus == 52, cut$c_plus == 53,
  abs(cut$tau - 0.3297297297297297) < 1e-14,
  abs(cut$tau * cut$P_plus + (1 - cut$tau) * cut$P_minus - 0.05) < 1e-14
)
print(cut$summary)

# Draw a fresh independent cutoff; tau is the probability of c_plus.
sampled <- qweighted_bernoulli(
  0.05, weights = 1:20, lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE
)
stopifnot(sampled$kappa %in% c(52, 53),
          identical(sampled$quantile, sampled$kappa))
print(sampled$kappa)

# The inclusive two-sided event counts an atom at the center only once.
center <- qweighted_bernoulli(
  0.75, weights = c(1, 1), lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE, draw = FALSE
)
stopifnot(center$exact, center$c_minus == 0, center$c_plus == 1,
          center$tau == 0.5, center$P_minus == 0.5, center$P_plus == 1)
print(center$summary)

# One-sided randomization works with unequal success probabilities.
one <- qweighted_bernoulli(
  0.05, weights = c(1, 2, 4), prob = c(0.2, 0.5, 0.8),
  lower.tail = FALSE, randomized = TRUE, draw = FALSE
)
stopifnot(one$exact, one$c_minus == 6, one$c_plus == 7,
          abs(one$tau - 0.375) < 1e-14,
          abs(one$tau * one$P_plus + (1 - one$tau) * one$P_minus - 0.05) < 1e-14)
print(one$summary)

# Short vectors: exact, including unequal Bernoulli probabilities.
q <- qweighted_bernoulli(c(0.5, 0.9, 0.95), c(1, 2, 4),
                         prob = c(0.2, 0.5, 0.8))
stopifnot(identical(unname(q), c(4, 6, 7)))
print(q)
q <- qweighted_bernoulli(c(0.025, 0.5, 0.975), 1:20)
stopifnot(identical(unname(q), c(53, 105, 157)))
print(q)

# The default max_length is 25; no method argument is needed.
short <- qweighted_bernoulli(0.95, sqrt(1:25), details = TRUE)
longer <- qweighted_bernoulli(0.95, sqrt(1:26), details = TRUE)
stopifnot(short$summary$exact, !longer$summary$exact)
print(short$summary)
print(longer$summary)

# A longer unstructured two-sided request has a conservative fallback.
bounded <- qweighted_bernoulli(
  0.05, weights = sqrt(1:26), lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE, draw = FALSE
)
stopifnot(!bounded$exact, bounded$c_minus == bounded$c_plus,
          bounded$tau == 1, bounded$P_plus <= 0.05 + 1e-14,
          bounded$summary$probability_type == "upper_bound")
print(bounded$summary)

# Exact structure can be used above the length threshold.
dominant <- qweighted_bernoulli(
  0.05, c(60, rep(1, 39)), lower.tail = FALSE, details = TRUE
)
stopifnot(dominant$summary$exact, dominant$quantile == 83)
print(dominant$summary)
q <- qweighted_bernoulli(0.95, rep(1, 1000))
stopifnot(q == stats::qbinom(0.95, 1000, 0.5))
print(q)

# General longer inputs: all applicable completed candidates are compared.
fit <- qweighted_bernoulli(
  c(0.05, 0.01, 1e-6), sqrt(1:100),
  lower.tail = FALSE, details = TRUE
)
print(fit$summary)
tab <- subset(fit$bounds, request == 1 & available)
tab <- tab[order(tab$quantile), ]
print(tab[, c("method", "quantile", "selected")])
stopifnot(fit$quantile[1] == min(tab$quantile))
print(subset(fit$bounds, !available,
             select = c("request", "method", "reason")))

# Keep very small upper-tail probabilities on the logarithmic scale.
q <- qweighted_bernoulli(
  log(1e-30), sqrt(1:26), lower.tail = FALSE, log.p = TRUE
)
stopifnot(is.finite(q))
print(q)

# Signed weights reduce to nonnegative weights plus a deterministic offset.
signed_w <- c(-3, 1, 2)
success <- c(0.2, 0.5, 0.8)
negative <- signed_w < 0
offset <- sum(signed_w[negative])
success[negative] <- 1 - success[negative]
q_signed <- offset + qweighted_bernoulli(
  0.95, abs(signed_w), prob = success
)
stopifnot(q_signed == 3)
print(q_signed)

dir.create("results", showWarnings = FALSE)
writeLines(c("All documented usage checks passed.", R.version.string,
             paste("Platform:", R.version$platform)),
           "results/usage_validation.txt")
cat("All documented usage checks passed.\n")
