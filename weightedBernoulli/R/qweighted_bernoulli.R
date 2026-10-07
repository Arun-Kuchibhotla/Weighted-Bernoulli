# Quantiles of independent weighted Bernoulli sums: ONE public function.
# X = sum(weights[i] * Bernoulli(prob[i])); nonnegative finite weights.
# Only base R and stats are needed. All implementation helpers are local.
#
# qweighted_bernoulli(p, weights, prob=0.5, max_length=25,
#                    lower.tail=TRUE, log.p=FALSE, details=FALSE,
#                    control=list())
#
# The ordinary p-quantile is inf{x: Pr(X <= x) >= p}. At p=0 and p=1
# return the support endpoints. lower.tail=FALSE supplies Pr(X > q);
# log.p=TRUE supplies its logarithm. Names and vector probability inputs work.
#
# Automatic dispatch (no method selector):
# 1. Remove zero and deterministic terms. Recognize binomial, fair rank,
#    affordable exact lattice, and repeated (weight, probability) count laws.
# 2. For at most max_length active terms, compute an exact finite-distribution
#    quantile by meet-in-the-middle if necessary. The default max_length is 25.
# 3. Otherwise evaluate every applicable candidate within its documented work
#    limits and return their minimum: concentration, upward grid, grouped
#    count envelopes, conditional moments, and conditional MGF-grid bounds.
# The count bounds currently require one common active success probability.
#
# A conservative result satisfies Pr(X > q) <= the requested upper-tail level.
# This is a STRICT tail. A returned bound need not be a support point.
# Exact means an exact finite-distribution algorithm in double precision,
# not symbolic or certified interval arithmetic. Numerically indistinguishable
# subset sums and CDF jumps retain the usual floating-point qualification.
#
# details=TRUE returns quantile, summary, bounds, and diagnostics. For a bound
# calculation, bounds records every attempted candidate, its availability,
# reason, partition, and selected flag; summary always selects its minimum.
# Exact requests have exact=TRUE and bypass unnecessary conservative work.
# The minimum is over the implemented finite candidate/partition/lambda set,
# not over all possible partitions or all probability inequalities.
#
# control defaults (integer limits; zero disables the indicated optional route):
# max_states=2^18: lattice/grid states and one MITM half-array.
# max_dp_work=2e7: cumulative lattice/grid updates.
# max_mitm_work=2e7: visited left-half entries in exact probability queries.
# max_iterations=80: conservative and structured-exact inversions.
# chernoff_iterations=60: optimized MGF-Chernoff evaluations per level.
# max_p2_B=10000, max_p2_work=200000: optional binomial P2 budgets.
# max_group_states=65536: exact-block/grouped pivot remainder states.
# max_count_states=100000: conditional joint count states.
# max_mgf_work=2e7: conditional MGF preparation/array budget per partition.
# mgf_points=37: finite grid, zero plus 36 log-spaced positive lambdas.
# Optional P2, grouped, count, and MGF budgets may be zero. For <=25 active
# terms, the promised exact branch floors MITM storage, bypasses the query-work cap,
# and allows 2200 support-search iterations; restrictive controls cannot silently
# switch it to an approximation. Larger user-selected max_length values are
# subject to the stated exact work limits and cause an error if unresolved.
#
# Examples:
# qweighted_bernoulli(c(0.025, 0.5, 0.975), weights=sqrt(1:25))
# qweighted_bernoulli(0.95, c(100, rep(1, 999)), details=TRUE)
# qweighted_bernoulli(0.05, seq(0.99, 1.01, length.out=1000),
#                    lower.tail=FALSE, details=TRUE)
# qweighted_bernoulli(log(1e-30), sqrt(1:100),
#                    lower.tail=FALSE, log.p=TRUE)
#
# Primary references: Bentkus-Dzindzalieta, doi:10.3150/14-BEJ603;
# Pinelis, doi:10.1214/EJP.v17-2026; Bentkus, arXiv:math/0410159;
# Keller-Klein, arXiv:2006.16834; Montgomery-Smith,
# doi:10.1090/S0002-9939-1990-1013975-0; Bardenet-Maillard,
# doi:10.3150/14-BEJ605. The accompanying note derives the count refinements.

qweighted_bernoulli <- function(
    p, weights, prob=0.5, max_length=25,
    lower.tail=TRUE, log.p=FALSE, details=FALSE, control=list()) {

  log_complement <- function(x) {
    answer <- numeric(length(x))
    small <- x < -log(2)
    answer[small] <- log1p(-exp(x[small]))
    answer[!small] <- log(-expm1(x[!small]))
    answer
  }

  probability_upper <- function(x) {
    ifelse(is.infinite(x) & x < 0, 0, pmin(1, pmax(2^-1074, exp(x))))
  }


  # Internal exact engines for the stand-alone weighted Bernoulli quantile.
  # These local helpers take positive weights as inputs;
  # probabilities are strictly between zero and one. The public wrapper handles
  # deterministic terms and endpoints. Each engine returns a data.frame with
  # quantile, exact, and method, or NULL when its resource limits do not permit
  # computation. An NA quantile is an unresolved request that needs a fallback.
  # lower_prob / upper_prob preserve unlogged input levels at CDF jumps.

  exact_target_ok <- function(cdf, survival, j, log_lower, log_upper,
                              lower_prob = NULL, upper_prob = NULL) {
    if (log_lower[j] <= -log(2)) {
      if (!is.null(lower_prob)) return(cdf >= lower_prob[j])
      return(log(cdf) >= log_lower[j])
    }
    if (!is.null(upper_prob)) return(survival <= upper_prob[j])
    log(survival) <= log_upper[j]
  }

  exact_result <- function(q, method) {
    data.frame(quantile = q, exact = TRUE, method = method,
               stringsAsFactors = FALSE)
  }

  # No approximate rounding of arbitrary scores is used to recognize a lattice.
  exact_lattice <- function(w) {
    integer_gcd <- function(a, b) {
      while (b != 0) {
        r <- a %% b
        a <- b
        b <- r
      }
      a
    }
    step <- min(w)
    z <- w / step
    if (all(is.finite(z)) && all(z == floor(z)) &&
        max(z) <= 2^52 && all(step * z == w) && sum(z) <= 2^52) {
      return(list(weights = z, step = step))
    }
    # Power-of-two scaling also recognizes, for example, weights (1, 1.5).
    # Try binary units without changing the represented input weights.
    exponent <- max(-1022, min(1023, floor(log2(max(w)))))
    step <- 2^exponent
    z <- w / step
    for (iteration in 0:52) {
      if (any(!is.finite(z)) || max(z) > 2^52) break
      if (all(z == floor(z)) && all(z > 0)) {
        g <- Reduce(integer_gcd, z)
        z <- z / g
        step <- step * g
        if (sum(z) <= 2^52 && all(z * step == w)) {
          return(list(weights = z, step = step))
        }
        break
      }
      z <- z * 2
      step <- step / 2
      if (step == 0) break
    }
    NULL
  }

  # Full lattice DP is also exposed for the wrapper's upward-rounded grid.
  # Its returned quantiles are in integer lattice coordinates.
  dp_quantiles <- function(w_integer, prob, log_lower, log_upper, control,
                           lower_prob = NULL, upper_prob = NULL) {
    ordering <- order(w_integer)
    w_integer <- w_integer[ordering]
    prob <- prob[ordering]
    total <- sum(w_integer)
    if (!is.finite(total) || total + 1 > control$max_states ||
        sum(1 + cumsum(w_integer)) > control$max_dp_work) return(NULL)
    if (any(w_integer < 1 | w_integer != floor(w_integer))) return(NULL)
    mass <- 1
    for (j in seq_along(w_integer)) {
      a <- w_integer[j]
      old <- seq_along(mass)
      next_mass <- numeric(length(mass) + a)
      next_mass[old] <- mass * (1 - prob[j])
      next_mass[old + a] <- next_mass[old + a] + mass * prob[j]
      mass <- next_mass
    }
    cdf <- cumsum(mass)
    # survival[k + 1] = P(sum > k), evaluated without 1 - CDF.
    survival <- c(rev(cumsum(rev(mass[-1L]))), 0)
    answer <- rep(NA_real_, length(log_lower))
    for (j in seq_along(answer)) {
      # Ordinary arithmetic is inappropriate for targets near probability
      # underflow. The public function may use a bound for these requests.
      if (min(log_lower[j], log_upper[j]) < -650) next
      ok <- function(k) exact_target_ok(cdf[k], survival[k], j,
        log_lower, log_upper, lower_prob, upper_prob)
      lo <- 1L
      hi <- length(mass)
      if (!ok(hi)) next
      while (lo < hi) {
        mid <- floor(lo + (hi - lo) / 2)
        if (ok(mid)) hi <- mid else lo <- mid + 1L
      }
      if (ok(lo) && (lo == 1L || !ok(lo - 1L))) answer[j] <- lo - 1
    }
    exact_result(answer, "lattice_dp")
  }

  # Correct the quantile routine's deliberate fuzz by checking the CDF jump.
  # The binary correction costs O(log(number of support points)) CDF calls.
  exact_standard_quantiles <- function(n, step, success, kind,
                                       log_lower, log_upper, control,
                                       lower_prob = NULL, upper_prob = NULL) {
    total <- if (kind == "binomial") n else n * (n + 1) / 2
    distribution <- if (kind == "binomial") {
      function(k, lower.tail, log.p) stats::pbinom(k, n, success,
        lower.tail = lower.tail, log.p = log.p)
    } else {
      function(k, lower.tail, log.p) stats::psignrank(k, n,
        lower.tail = lower.tail, log.p = log.p)
    }
    quantile <- if (kind == "binomial") {
      function(p, lower.tail) stats::qbinom(p, n, success,
        lower.tail = lower.tail, log.p = TRUE)
    } else {
      function(p, lower.tail) stats::qsignrank(p, n,
        lower.tail = lower.tail, log.p = TRUE)
    }
    if (kind == "binomial" && success == 0.5 && n <= 52L) {
      row <- 1
      for (j in seq_len(n)) row <- c(row, 0) + c(0, row)
      cdf <- cumsum(row) * 2^-n
      survival <- c(rev(cumsum(rev(row[-1L]))), 0) * 2^-n
      distribution <- function(k, lower.tail, log.p) {
        at <- pmin(n, pmax(0, floor(k))) + 1L
        value <- if (lower.tail) cdf[at] else survival[at]
        value[k < 0] <- if (lower.tail) 0 else 1
        value[k >= n] <- if (lower.tail) 1 else 0
        if (log.p) log(value) else value
      }
    }
    answer <- rep(NA_real_, length(log_lower))
    for (j in seq_along(answer)) {
      lower <- log_lower[j] <= -log(2)
      ok <- if (lower) {
        if (is.null(lower_prob)) {
          function(k) distribution(k, TRUE, TRUE) >= log_lower[j]
        } else {
          function(k) distribution(k, TRUE, FALSE) >= lower_prob[j]
        }
      } else {
        if (is.null(upper_prob)) {
          function(k) distribution(k, FALSE, TRUE) <= log_upper[j]
        } else {
          function(k) distribution(k, FALSE, FALSE) <= upper_prob[j]
        }
      }
      q <- quantile(if (lower) log_lower[j] else log_upper[j], lower)
      if (!is.finite(q)) next
      q <- min(total, max(0, q))
      lo <- 0
      hi <- total
      if (!ok(q)) {
        lo <- q + 1
      } else if (q > 0 && ok(q - 1)) {
        hi <- q - 1
      } else {
        lo <- hi <- q
      }
      if (lo > hi) next
      for (iteration in seq_len(control$max_iterations)) {
        if (lo >= hi) break
        mid <- floor(lo + (hi - lo) / 2)
        if (ok(mid)) hi <- mid else lo <- mid + 1
      }
      if (lo == hi && ok(lo) && (lo == 0 || !ok(lo - 1))) {
        answer[j] <- step * lo
      }
    }
    exact_result(answer, kind)
  }

  exact_half_distribution <- function(w, prob) {
    values <- 0
    mass <- 1
    log_mass <- 0
    for (j in seq_along(w)) {
      values <- c(values, values + w[j])
      mass <- c(mass * (1 - prob[j]), mass * prob[j])
      log_mass <- c(log_mass + log1p(-prob[j]), log_mass + log(prob[j]))
    }
    ordering <- order(values)
    list(values = values[ordering], mass = mass[ordering],
         log_mass = log_mass[ordering])
  }

  exact_mitm_quantiles <- function(w, prob, log_lower, log_upper, control,
                                  lower_prob = NULL, upper_prob = NULL) {
    n <- length(w)
    if (n > min(52, control$mitm_B) ||
        2^ceiling(n / 2) > control$max_states) return(NULL)
    m <- floor(n / 2)
    left <- exact_half_distribution(w[seq_len(m)], prob[seq_len(m)])
    right_indices <- seq.int(m + 1L, n)
    right <- exact_half_distribution(w[right_indices], prob[right_indices])
    a <- left$values
    b <- right$values
    nb <- length(b)
    prefix <- c(0, cumsum(right$mass))
    suffix <- c(rev(cumsum(rev(right$mass))), 0)
    use_log <- any(pmin(log_lower, log_upper) < -650)
    log_add <- function(x, y) {
      z <- max(x, y)
      if (!is.finite(z)) return(z)
      z + log1p(exp(-abs(x-y)))
    }
    log_sum <- function(x) {
      z <- max(x)
      if (!is.finite(z)) return(z)
      z + log(sum(exp(x-z)))
    }
    if (use_log) {
      log_prefix <- log_suffix <- rep(-Inf, nb + 1L)
      for (j in seq_len(nb))
        log_prefix[j + 1L] <- log_add(log_prefix[j], right$log_mass[j])
      for (j in rev(seq_len(nb)))
        log_suffix[j] <- log_add(log_suffix[j + 1L], right$log_mass[j])
    }
    used <- 0
    maximum <- max(a) + max(b)
    cache <- new.env(parent = emptyenv())

    pairs <- function(x, strict = FALSE) {
      used <<- used + length(a)
      if (used > control$max_mitm_work) return(NULL)
      j <- findInterval(x - a, b)
      bad <- logical(length(a))
      ii <- which(j > 0L)
      v <- a[ii] + b[j[ii]]
      bad[ii] <- if (strict) v >= x else v > x
      ii <- which(j < nb)
      v <- a[ii] + b[j[ii] + 1L]
      bad[ii] <- bad[ii] | if (strict) v < x else v <= x
      ii <- which(bad)
      if (length(ii)) {
        # Repair subtraction/addition rounding disagreements. This also
        # moves to the correct edge of a repeated subset-sum value.
        lower <- integer(length(ii))
        upper <- rep.int(nb + 1L, length(ii))
        active <- seq_along(ii)
        while (length(active)) {
          mid <- floor((lower[active] + upper[active]) / 2)
          v <- a[ii[active]] + b[mid]
          good <- if (strict) v < x else v <= x
          lower[active[good]] <- mid[good]
          upper[active[!good]] <- mid[!good]
          active <- which(upper - lower > 1L)
        }
        j[ii] <- lower
      }
      list(cdf = sum(left$mass * prefix[j + 1L]),
           survival = sum(left$mass * suffix[j + 1L]),
           log_cdf = if (use_log)
             min(0, log_sum(left$log_mass + log_prefix[j + 1L])) else NULL,
           log_survival = if (use_log)
             min(0, log_sum(left$log_mass + log_suffix[j + 1L])) else NULL,
           j = j)
    }

    answer <- rep(NA_real_, length(log_lower))
    for (k in seq_along(answer)) {
      lower <- log_lower[k] <= -log(2)
      target <- if (lower) {
        if (is.null(lower_prob)) log_lower[k] else lower_prob[k]
      } else {
        if (is.null(upper_prob)) log_upper[k] else upper_prob[k]
      }
      key <- paste0(if (lower) "L" else "U", sprintf("%.17g", target))
      if (exists(key, envir = cache, inherits = FALSE)) {
        answer[k] <- get(key, envir = cache, inherits = FALSE)
        next
      }
      ok <- function(z) {
        if (use_log && min(log_lower[k], log_upper[k]) < -650) {
          if (lower) return(z$log_cdf >= log_lower[k])
          return(z$log_survival <= log_upper[k])
        }
        exact_target_ok(z$cdf, z$survival, k,
          log_lower, log_upper, lower_prob, upper_prob)
      }
      low_query <- pairs(0)
      if (is.null(low_query)) break
      if (ok(low_query)) {
        lo <- hi <- 0
      } else {
        lo <- 0
        hi <- maximum
        for (iteration in seq_len(control$max_iterations)) {
          if (lo == hi) break
          mid <- lo + (hi - lo) / 2
          if (mid == hi) mid <- lo
          z <- pairs(mid)
          if (is.null(z)) break
          if (ok(z)) {
            ii <- which(z$j > 0L)
            if (!length(ii)) break
            hi <- max(a[ii] + b[z$j[ii]])
          } else {
            ii <- which(z$j < nb)
            if (!length(ii)) break
            lo <- min(a[ii] + b[z$j[ii] + 1L])
          }
          if (lo > hi) break
        }
      }
      if (lo != hi) next
      plus <- pairs(lo)
      minus <- pairs(lo, strict = TRUE)
      if (is.null(plus) || is.null(minus)) break
      if (ok(plus) && !ok(minus)) {
        answer[k] <- lo
        assign(key, lo, envir = cache)
      }
    }
    exact_result(answer, "meet_in_the_middle")
  }

  # Exact compression of repeated (weight, probability) pairs into independent
  # binomial counts. Integrate the largest group by its CDF and enumerate only
  # the other count states. A NULL result means that the structure/work limits
  # do not permit this route; unresolved numerical queries return NA.
  structured_count_quantiles <- function(
      w, prob, log_lower, log_upper, control,
      lower_prob = NULL, upper_prob = NULL) {
    n <- length(w)
    if (n < 2L) return(NULL)
    key <- paste(sprintf("%.17g", w), sprintf("%.17g", prob), sep = ":")
    group <- match(key, unique(key))
    ng <- max(group)
    m <- tabulate(group, nbins = ng)
    if (all(m == 1L)) return(NULL)
    first <- match(seq_len(ng), group)
    gw <- w[first]
    gp <- prob[first]
    pivot <- which.max(m)
    log_states <- sum(log1p(m[-pivot]))
    if (log_states > log(control$max_group_states) + 1e-12) return(NULL)
    states <- prod(m[-pivot] + 1)
    if (!is.finite(states) || states > control$max_group_states) return(NULL)

    log_sum <- function(x) {
      a <- max(x)
      if (!is.finite(a)) return(a)
      a + log(sum(exp(x - a)))
    }
    count_mass <- function(size, success) {
      k <- seq.int(0, size)
      if (success == 0.5 && size <= 52L) {
        a <- 1
        for (i in seq_len(size)) a <- c(a, 0) + c(0, a)
        mass <- a * 2^-size
        return(list(mass = mass, log_mass = log(mass)))
      }
      list(mass = stats::dbinom(k, size, success),
           log_mass = stats::dbinom(k, size, success, log = TRUE))
    }
    values <- 0
    mass <- 1
    log_mass <- 0
    for (j in setdiff(seq_len(ng), pivot)) {
      counts <- seq.int(0, m[j])
      tab <- count_mass(m[j], gp[j])
      old <- length(values)
      values <- as.vector(outer(values, gw[j] * counts, "+"))
      mass <- rep(mass, times = length(counts)) * rep(tab$mass, each = old)
      log_mass <- rep(log_mass, times = length(counts)) +
        rep(tab$log_mass, each = old)
    }
    pv <- gw[pivot] * seq.int(0, m[pivot])
    nv <- length(pv)
    maximum <- max(values) + pv[nv]
    if (!is.finite(maximum)) return(NULL)
    if (gp[pivot] == 0.5 && m[pivot] <= 52L) {
      pm <- count_mass(m[pivot], gp[pivot])$mass
      prefix <- c(0, cumsum(pm))
      suffix <- c(rev(cumsum(rev(pm))), 0)
      log_prefix <- log(prefix)
      log_suffix <- log(suffix)
    } else {
      thresholds <- seq.int(-1, m[pivot])
      prefix <- stats::pbinom(thresholds, m[pivot], gp[pivot])
      suffix <- stats::pbinom(thresholds, m[pivot], gp[pivot], lower.tail = FALSE)
      log_prefix <- stats::pbinom(thresholds, m[pivot], gp[pivot], log.p = TRUE)
      log_suffix <- stats::pbinom(thresholds, m[pivot], gp[pivot],
                                 lower.tail = FALSE, log.p = TRUE)
    }

    pair_query <- function(x, strict = FALSE) {
      j <- findInterval(x - values, pv)
      bad <- logical(length(values))
      ii <- which(j > 0L)
      v <- values[ii] + pv[j[ii]]
      bad[ii] <- if (strict) v >= x else v > x
      ii <- which(j < nv)
      v <- values[ii] + pv[j[ii] + 1L]
      bad[ii] <- bad[ii] | if (strict) v < x else v <= x
      ii <- which(bad)
      if (length(ii)) {
        lo <- integer(length(ii))
        hi <- rep.int(nv + 1L, length(ii))
        active <- seq_along(ii)
        while (length(active)) {
          mid <- floor((lo[active] + hi[active]) / 2)
          v <- values[ii[active]] + pv[mid]
          good <- if (strict) v < x else v <= x
          lo[active[good]] <- mid[good]
          hi[active[!good]] <- mid[!good]
          active <- which(hi - lo > 1L)
        }
        j[ii] <- lo
      }
      list(cdf = min(1, sum(mass * prefix[j + 1L])),
           survival = min(1, sum(mass * suffix[j + 1L])),
           log_cdf = min(0, log_sum(log_mass + log_prefix[j + 1L])),
           log_survival = min(0, log_sum(log_mass + log_suffix[j + 1L])),
           j = j)
    }

    target_ok <- function(z, i) {
      if (log_lower[i] <= -log(2)) {
        if (!is.null(lower_prob) && log_lower[i] >= -650)
          return(z$cdf >= lower_prob[i])
        return(z$log_cdf >= log_lower[i])
      }
      if (!is.null(upper_prob) && log_upper[i] >= -650)
        return(z$survival <= upper_prob[i])
      z$log_survival <= log_upper[i]
    }
    ans <- rep(NA_real_, length(log_lower))
    for (i in seq_along(ans)) {
      lo <- 0
      hi <- maximum
      if (target_ok(pair_query(0), i)) hi <- 0
      for (iteration in seq_len(control$max_iterations)) {
        if (lo == hi) break
        mid <- lo + (hi - lo) / 2
        if (mid == hi) mid <- lo
        z <- pair_query(mid)
        if (target_ok(z, i)) {
          ii <- which(z$j > 0L)
          if (!length(ii)) break
          hi <- max(values[ii] + pv[z$j[ii]])
        } else {
          ii <- which(z$j < nv)
          if (!length(ii)) break
          lo <- min(values[ii] + pv[z$j[ii] + 1L])
        }
        if (lo > hi) break
      }
      if (lo == hi && target_ok(pair_query(lo), i) &&
          !target_ok(pair_query(lo, strict = TRUE), i)) ans[i] <- lo
    }
    exact_result(ans, "binomial_blocks")
  }

  # Nested inside the public weighted-Bernoulli quantile function.
  # Inputs: positive normalized w; 0 < prob[i] < 1; vector log_beta <= 0.
  # Returns one list per log_beta with quantile, method, log_tail_bound,
  # tail_bound, and a data.frame of the individual candidates.
  # Every finite candidate bounds the STRICT tail P(sum(w*Bernoulli(prob)) > q).
  bound_candidates <- function(w, prob, log_beta, control) {
    eps <- .Machine$double.eps
    tiny <- 2^-1074
    n <- length(w)
    total <- sum(w)
    center <- sum(w*prob)
    variance <- sum((w*sqrt(prob)*sqrt(1-prob))^2)*(1+128*eps) + n*tiny
    sd <- sqrt(variance)
    width2 <- sum(w*w)*(1+128*eps) + n*tiny
    y <- max(w*(1-prob))*(1+64*eps) + tiny
    symmetric <- all(prob == 0.5)
    margin <- max(256*eps*total, 8*tiny)
    log_p <- log(prob)
    log_q <- log1p(-prob)
    log_top <- sum(log_p)
    p2_remaining <- control$max_p2_work
    p2_prepared <- NULL
    p2_cache <- new.env(parent=emptyenv())
    k_prepared <- NULL
    log_C_BD <- -log(4)-stats::pnorm(sqrt(2), lower.tail=FALSE, log.p=TRUE)
    log_C_P <- log(5*sqrt(2*base::pi*exp(1))*(2*stats::pnorm(1)-1))

    log_add <- function(a,b) pmax(a,b)+log1p(exp(-abs(a-b)))
    log1pexp <- function(a) pmax(a,0)+log1p(exp(-abs(a)))
    log_Q <- function(x) {
      log_add(stats::pnorm(x, lower.tail=FALSE, log.p=TRUE),
              log_C_P+stats::dnorm(x, log=TRUE)-log(9+x*x))
    }
    inverse_Q <- function(target) {
      lo <- 0
      hi <- max(1, stats::qnorm(target-log_C_BD, lower.tail=FALSE, log.p=TRUE))
      for (j in seq_len(8L)) {
        if (log_Q(hi) <= target) break
        hi <- 2*hi
      }
      if (!is.finite(hi) || log_Q(hi) > target) return(Inf)
      for (j in seq_len(min(60,control$max_iterations))) {
        mid <- lo+(hi-lo)/2
        if (mid==lo || mid==hi) break
        if (log_Q(mid)>target) lo <- mid else hi <- mid
      }
      hi
    }

    # K(a,x; l1,l2) = inf_b {sum(abs(b)) + x*sqrt(sum((a-b)^2))}.
    # Montgomery-Smith's deterministic-part argument also applies to any
    # Gaussian Rademacher comparison. The minimizing residual caps the largest
    # coefficients at r: r^2=sum(a[(k+1):n]^2)/(x^2-k), for decreasing a.
    # w is already increasing. Cache its reverse cumulative second moments;
    # only k < x^2 cells can be relevant. The returned head/norm are evaluated
    # from an actual decomposition, so an inexact cap affects tightness only.
    gaussian_k_split <- function(x) {
      if (!is.finite(x) || x<=0 || x>=sqrt(n) ||
          x*max(w)<=sqrt(width2)) return(NULL)
      if (is.null(k_prepared))
        k_prepared <<- list(a=rev(w), tail2=rev(cumsum(w*w)))
      t2 <- x*x
      ii <- seq_len(min(n-1,ceiling(t2)-1)+1)
      cap <- sqrt(k_prepared$tail2[ii]/(t2-(ii-1)))
      jj <- which(is.finite(cap) & cap>=k_prepared$a[ii])
      if (!length(jj)) return(NULL)
      cap <- cap[jj[1L]]
      if (cap<=0 || cap>=max(w)) return(NULL)
      residual <- pmin(w,cap)
      head <- (sum(pmax(w-cap,0))/2)*(1+128*eps)+n*tiny
      residual_variance <- sum((residual/2)^2)*(1+128*eps)+n*tiny
      residual_sd <- sqrt(residual_variance)*(1+16*eps)+tiny
      radius <- head+x*residual_sd
      # Suppress numerical ties with the original Gaussian comparison.
      if (!is.finite(radius) || radius>=sd*x-128*eps*total) return(NULL)
      list(head=head,sd=residual_sd,radius=radius)
    }

    # Bentkus (2004), Lemmas 4.4--4.5, https://arxiv.org/pdf/math/0410159.
    # X_i=w_i*(Bernoulli(prob_i)-prob_i) <= y, with total variance <= variance.
    # The centered comparison variable is step*(Bin(n,pc)-n*pc).
    prepare_p2 <- function() {
      ratio <- variance/(n*y*y)
      pc <- (ratio/(1+ratio))*(1+16*eps)+tiny
      if (!is.finite(pc) || pc<=0 || pc>=1) return(NULL)
      step <- y/(1-pc)*(1+16*eps)
      if (!is.finite(step)) return(NULL)
      k <- seq.int(n,0)
      pmf <- stats::dbinom(k,n,pc,log=TRUE)
      tail_logs <- deltas <- variances <- numeric(n+1L)
      tail_log <- pmf[1L]
      delta <- v <- 0
      for (i in seq_along(k)) {
        if (i>1L) {
          combined <- log_add(pmf[i],tail_log)
          wt <- exp(tail_log-combined)
          wa <- exp(pmf[i]-combined)
          old_mean <- delta+1
          v <- wt*v+wt*wa*old_mean^2
          delta <- wt*old_mean
          tail_log <- combined
        }
        tail_logs[i] <- tail_log
        deltas[i] <- delta
        variances[i] <- v
      }
      list(p=pc,step=step,k=k,tail_log=tail_logs,delta=deltas,
           variance=variances,guard=256*eps*(n+1))
    }

    # Analytic optimization within each binomial support cell. Negative shifts
    # are included in the k=0 cell; preparation is shared across target levels.
    inverse_p2 <- function(target) {
      if (n>control$max_p2_B || control$max_p2_B==0) return(NULL)
      key <- sprintf("%.17g",target)
      if (exists(key,envir=p2_cache,inherits=FALSE))
        return(get(key,envir=p2_cache,inherits=FALSE))
      if (is.null(p2_prepared)) {
        if (2*(n+1)>p2_remaining) return(NULL)
        p2_remaining <<- p2_remaining-(n+1)
        p2_prepared <<- prepare_p2()
        if (is.null(p2_prepared)) return(NULL)
      }
      if (n+1>p2_remaining) return(NULL)
      p2_remaining <<- p2_remaining-(n+1)
      z <- p2_prepared
      if (target<=n*log(z$p)+z$guard) return(NULL)
      k <- z$k
      d <- z$delta
      v <- z$variance
      lr <- z$tail_log+z$guard-target
      offset <- rep(-1,length(k))
      offset[k==0] <- NA_real_
      ii <- which(lr>0 & v>0)
      if (length(ii)) {
        a <- lr[ii]
        ld <- a+log(-expm1(-a))
        lg <- (log(v[ii])-ld)/2
        interior <- k[ii]==0 | lg<log1p(d[ii])
        jj <- ii[interior]
        offset[jj] <- pmin(0,d[jj]-exp(lg[interior]))
      }
      moment <- (d-offset)^2+v
      lm <- z$tail_log+log(moment)+z$guard
      lradius <- (lm-target)/2
      good <- is.finite(moment) & moment>0 & is.finite(lradius) &
        lradius<log(.Machine$double.xmax)
      good[is.na(good)] <- FALSE
      ii <- which(good)
      if (!length(ii)) return(NULL)
      qb <- k[ii]+offset[ii]+exp(lradius[ii])
      j <- ii[which.min(qb)]
      ans <- list(radius=z$step*(min(qb)-n*z$p),
                  shift=z$step*(k[j]+offset[j]-n*z$p),
                  log_moment=2*log(z$step)+lm[j])
      assign(key,ans,envir=p2_cache)
      ans
    }

    # Optimize the Chernoff quantile directly: inf_{lambda>0}
    # {log E exp(lambda*G)-log(beta)}/lambda. Every evaluated lambda is valid.
    # The entropy equation lambda*K'(lambda)-K(lambda)=-log(beta) is monotone.
    chernoff <- function(target) {
      lambda <- min(1/max(w),sqrt(-2*target/variance))
      lo <- 0
      hi <- Inf
      best <- total
      certificate <- NULL
      for (j in seq_len(control$chernoff_iterations)) {
        z <- lambda*w
        if (any(!is.finite(z))) break
        g <- log_add(log_q,log_p+z)
        small <- z<1
        g[small] <- log1p(prob[small]*expm1(z[small]))
        K <- sum(g)
        k <- log_add(log_p,log_q-z)
        entropy_terms <- -k-z*stats::plogis(log_q-log_p-z)
        entropy_terms[small] <- z[small]*
          stats::plogis(log_p[small]-log_q[small]+z[small])-g[small]
        entropy <- sum(entropy_terms)
        guard <- 128*eps*(n+abs(K)+abs(target)+1)
        q <- (K-target+2*guard)/lambda
        if (is.finite(q) && q<best) {
          best <- q
          certificate <- list(lambda=lambda,K=K,guard=guard)
        }
        if (!is.finite(entropy)) break
        if (entropy < -target) lo <- lambda else hi <- lambda
        if (is.finite(hi) && hi-lo<=sqrt(eps)*hi) break
        next_lambda <- if (is.finite(hi)) lo+(hi-lo)/2 else 2*lambda
        if (!is.finite(next_lambda) || next_lambda==lambda) break
        lambda <- next_lambda
      }
      list(quantile=best,certificate=certificate)
    }

    lapply(log_beta,function(lb) {
      row <- function(method,q,log_tail) {
        data.frame(method=method,quantile=q,log_tail_bound=log_tail,
                   stringsAsFactors=FALSE)
      }
      # The upward allowance also covers the final sum of represented weights.
      # Without it, a rounded-down support endpoint can falsely exclude its atom.
      support <- total+margin
      candidates <- list(row("support",support,-Inf))
      if (lb==0) candidates[[2L]] <- row("support",0,0)
      if (is.finite(lb) && lb<0) {
        target <- lb-64*eps*(1+abs(lb))
        unavailable <- function(name) {
          candidates[[length(candidates)+1L]] <<- row(name,NA_real_,NA_real_)
        }
        add <- function(name,raw,log_tail) {
          if (!is.finite(raw)) return(unavailable(name))
          q <- max(0,raw+margin)
          if (q>=total) {
            candidates[[length(candidates)+1L]] <<- row(name,support,-Inf)
            return(invisible(NULL))
          }
          # Recheck at a slightly smaller threshold after the rounding allowance.
          certificate <- log_tail(max(0,q-margin/2))
          if (!is.finite(certificate) || certificate>lb)
            return(unavailable(name))
          candidates[[length(candidates)+1L]] <<- row(name,q,certificate)
          invisible(NULL)
        }
        atom_guard <- 64*eps*(1+abs(log_top))
        if (log_top+atom_guard<=lb)
          add("endpoint",total-min(w),function(q) log_top+atom_guard)

        z <- chernoff(target)
        if (is.null(z$certificate)) unavailable("chernoff") else {
          cert <- z$certificate
          add("chernoff",z$quantile,function(q)
            cert$K-cert$lambda*q+cert$guard)
        }
        add("hoeffding",center+sqrt(-width2*target/2),function(q)
          if (q>center) -2*((q-center)/sqrt(width2))^2 else 0)
        log_cantelli_radius <- (log(variance)-target+log(-expm1(target)))/2
        add("cantelli",center+exp(log_cantelli_radius),function(q)
          if (q>center) -log1pexp(2*log(q-center)-log(variance)) else 0)
        t <- -y*target/3
        add("bernstein",center+t+sqrt(t*t-2*variance*target),function(q) {
          r <- q-center
          if (r>0) -r*r/(2*(variance+y*r/3)) else 0
        })
        if (symmetric) {
          # Bentkus--Dzindzalieta, Theorem 1.1: https://arxiv.org/pdf/1307.3451.
          x <- stats::qnorm(target-log_C_BD,lower.tail=FALSE,log.p=TRUE)
          add("bentkus_gaussian",center+sd*x,function(q)
            log_C_BD+stats::pnorm((q-center)/sd,lower.tail=FALSE,log.p=TRUE))
          k_bd <- gaussian_k_split(x)
          if (is.null(k_bd)) unavailable("bentkus_gaussian_k") else
            add("bentkus_gaussian_k",center+k_bd$radius,function(q)
              log_C_BD+stats::pnorm((q-center-k_bd$head)/k_bd$sd,
                                    lower.tail=FALSE,log.p=TRUE))
          # Pinelis, Theorem 1.1, (1.7)--(1.8): https://arxiv.org/html/1007.2137v3.
          x <- inverse_Q(target)
          add("pinelis",center+sd*x,function(q) log_Q((q-center)/sd))
          k_pi <- gaussian_k_split(x)
          if (is.null(k_pi)) unavailable("pinelis_k") else
            add("pinelis_k",center+k_pi$radius,function(q)
              log_Q((q-center-k_pi$head)/k_pi$sd))
          add("symmetric_chebyshev",center+exp((log(variance)-log(2)-target)/2),
              function(q) if (q>center) log(variance)-log(2)-2*log(q-center) else 0)
          # Keller--Klein, Theorem 1.2: https://arxiv.org/pdf/2006.16834.
          # Symmetry halves their strict two-sided tail bound of 1/2.
          if (lb>=-log(4)) add("tomaszewski",center+sd,function(q)
            if (q-center>=sd) -log(4) else 0)
          if (lb>=-log(2)) add("symmetry",center,function(q)
            if (q>=center) -log(2) else 0)
        }
        p2 <- inverse_p2(target)
        if (is.null(p2)) unavailable("binomial_p2") else {
          add("binomial_p2",center+p2$radius,function(q) {
            gap <- q-center-p2$shift
            if (gap>0) p2$log_moment-2*log(gap) else 0
          })
        }
      }
      tab <- do.call(rbind,candidates)
      tab$available <- is.finite(tab$quantile) & !is.na(tab$log_tail_bound) &
        tab$log_tail_bound<=lb
      valid <- which(tab$available)
      winner <- valid[which.min(tab$quantile[valid])]
      tab$selected <- seq_len(nrow(tab))==winner
      log_bound <- tab$log_tail_bound[winner]
      list(quantile=tab$quantile[winner],method=tab$method[winner],
           log_tail_bound=log_bound,
           tail_bound=if (log_bound==-Inf) 0 else max(tiny,exp(log_bound)),
           candidates=tab)
    })
  }
  grid_quantiles <- function(w, prob, log_lower, log_upper, control,
                             lower_prob=NULL, upper_prob=NULL) {
    n <- length(w)
    states <- min(control$max_states, floor(control$max_dp_work/n))
    # Every positive weight costs at least one grid unit when rounded upward.
    if (states <= n+1) return(NULL)
    units <- states-1
    exponent <- ceiling(log2(sum(w))-log2(units-n))
    step <- 2^max(-1074, min(1023, exponent))
    repeat {
      integer_weights <- ceiling(w/step)
      if (any(!is.finite(integer_weights))) return(NULL)
      # Explicitly retain the direction of the coupling after multiplication.
      too_small <- integer_weights*step < w
      integer_weights[too_small] <- integer_weights[too_small]+1
      if (sum(integer_weights) <= units) break
      step <- 2*step
      if (!is.finite(step)) return(NULL)
    }
    proxy <- step*integer_weights
    if (any(!is.finite(proxy)) || any(proxy < w)) return(NULL)
    z <- dp_quantiles(integer_weights, prob, log_lower, log_upper, control,
                      lower_prob=lower_prob, upper_prob=upper_prob)
    if (is.null(z)) return(NULL)
    z$quantile <- step*z$quantile
    z$exact <- all(proxy == w)
    error <- sum(proxy-w)
    if (error > 0) error <- error*(1+64*.Machine$double.eps)+
      64*.Machine$double.eps*sum(w)
    z$error <- error
    z
  }

  grouped_quantiles <- function(
      p, weights, prob=0.5, groups=NULL, exact=integer(),
      lower.tail=TRUE, log.p=FALSE, details=FALSE, control=list()) {

    scalar_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
    if (!is.numeric(weights) || anyNA(weights) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be a finite nonnegative numeric vector")
    }
    if (!is.numeric(prob) || length(prob) != 1L || is.na(prob) ||
        !is.finite(prob) || prob < 0 || prob > 1) {
      stop("prob must be one common Bernoulli success probability in [0,1]")
    }
    if (!scalar_flag(lower.tail) || !scalar_flag(log.p) ||
        !scalar_flag(details)) {
      stop("lower.tail, log.p, and details must be single logical values")
    }
    if (!is.numeric(p) || anyNA(p) ||
        (!log.p && any(!is.finite(p) | p < 0 | p > 1)) ||
        (log.p && any(p > 0))) {
      stop("p must contain probabilities in [0,1], or log probabilities <= 0")
    }
    defaults <- list(max_states=65536, max_iterations=64)
    if (!is.list(control) ||
        (length(control) && (is.null(names(control)) ||
         any(!names(control) %in% names(defaults)) ||
         anyDuplicated(names(control))))) {
      stop("control must be a named list containing max_states or max_iterations")
    }
    for (nm in names(control)) defaults[[nm]] <- control[[nm]]
    for (nm in names(defaults)) {
      value <- defaults[[nm]]
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
          value < 1 || value != floor(value)) {
        stop(paste0("control$", nm, " must be a positive finite integer"))
      }
    }
    if (defaults$max_states > .Machine$integer.max) {
      stop("control$max_states must not exceed .Machine$integer.max")
    }
    control <- defaults
    total <- sum(weights)
    if (!is.finite(total)) stop("sum(weights) must be finite")

    check_indices <- function(x, what) {
      if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
          any(x != floor(x) | x < 1 | x > length(weights)) ||
          anyDuplicated(x)) {
        stop(paste0(what, " must contain distinct valid weight indices"))
      }
      as.integer(x)
    }
    exact <- check_indices(exact, "exact")
    if (!is.null(groups)) {
      if (!is.list(groups)) stop("groups must be a list of index vectors")
      if (length(exact)) stop("use either explicit groups or exact, not both")
      groups <- lapply(groups, check_indices, what="each group")
      used <- unlist(groups, use.names=FALSE)
      if (anyDuplicated(used)) stop("groups must be disjoint")
    }

    log_complement <- function(x) {
      ans <- numeric(length(x))
      small <- x < -log(2)
      ans[small] <- log1p(-exp(x[small]))
      ans[!small] <- log(-expm1(x[!small]))
      ans
    }
    log_input <- if (log.p) p else log(p)
    log_other <- log_complement(log_input)
    log_lower <- if (lower.tail) log_input else log_other
    log_upper <- if (lower.tail) log_other else log_input
    numeric_lower <- if (!log.p) if (lower.tail) p else 1-p else NULL
    numeric_upper <- if (!log.p) if (lower.tail) 1-p else p else NULL
    use_lower <- log_lower <= -log(2)

    active <- which(weights > 0)
    singleton_result <- function(value) {
      answer <- rep(value, length(p))
      if (!details) return(answer)
      list(quantile=answer, lower=answer, deterministic_error=0,
        groups=list(), group_summary=data.frame(), pivot=NA_integer_,
        enumerated_states=1, total_count_states=1, exact_grouping=TRUE,
        summary=data.frame(p=p, lower=answer, quantile=answer,
          bracket_width=rep(0, length(p)), search_error=rep(0, length(p)),
          resolved=rep(TRUE, length(p))))
    }
    if (!length(active) || prob == 0) return(singleton_result(0))
    if (prob == 1) return(singleton_result(total))

    gap <- function(indices) {
      a <- sort(weights[indices], decreasing=TRUE)
      h <- length(a) %/% 2L
      if (h == 0L || a[1L] == a[length(a)]) return(0)
      max(0, sum(a[seq_len(h)]) - sum(tail(a, h)))
    }
    log_states <- function(blocks) {
      m <- lengths(blocks)
      if (!length(m)) return(0)
      sum(log1p(m)) - log1p(max(m))
    }
    log_budget <- log(control$max_states)

    if (is.null(groups)) {
      kept <- setdiff(active, exact)
      kept <- kept[order(weights[kept], decreasing=TRUE)]
      groups <- if (length(kept)) list(kept) else list()
      exact <- intersect(exact, active)
      groups <- c(groups, lapply(exact, function(i) i))
      if (log_states(groups) > log_budget + 1e-12) {
        stop("the requested exact indices exceed control$max_states")
      }

      # Each candidate is a split between adjacent sorted coefficients.
      # The choice minimizes the resulting deterministic error greedily.
      repeat {
        m <- lengths(groups)
        all_log <- sum(log1p(m))
        best_gain <- 0
        best_cost <- Inf
        best_group <- best_split <- NA_integer_
        for (j in seq_along(groups)) {
          a <- weights[groups[[j]]]
          n <- length(a)
          if (n < 2L || a[1L] == a[n]) next
          k <- seq_len(n - 1L)
          other_max <- if (length(m) == 1L) 0 else max(m[-j])
          cost <- all_log - log1p(n) + log1p(k) + log1p(n-k) -
            log1p(pmax(other_max, k, n-k))
          feasible <- cost <= log_budget + 1e-12
          if (!any(feasible)) next
          z <- c(0, cumsum(a))
          h1 <- k %/% 2L
          h2 <- (n-k) %/% 2L
          left_gap <- z[h1+1L] - z[k+1L] + z[k-h1+1L]
          right_gap <- z[k+h2+1L] - z[k+1L] - z[n+1L] + z[n-h2+1L]
          gain <- gap(groups[[j]]) - pmax(0, left_gap) - pmax(0, right_gap)
          gain[!feasible] <- -Inf
          candidates <- which(gain == max(gain))
          at <- candidates[which.min(cost[candidates])]
          if (gain[at] > best_gain ||
              (gain[at] == best_gain && gain[at] > 0 && cost[at] < best_cost)) {
            best_gain <- gain[at]
            best_cost <- cost[at]
            best_group <- j
            best_split <- k[at]
          }
        }
        if (is.na(best_group)) break
        old <- groups[[best_group]]
        groups[[best_group]] <- old[seq_len(best_split)]
        groups <- append(groups, list(old[-seq_len(best_split)]), best_group)
      }
    } else {
      used <- unlist(groups, use.names=FALSE)
      groups <- lapply(groups, function(g) intersect(g, active))
      groups <- groups[lengths(groups) > 0L]
      groups <- c(groups, lapply(setdiff(active, used), function(i) i))
    }
    groups <- lapply(groups, function(g) g[order(weights[g], decreasing=TRUE)])
    m <- lengths(groups)
    pivot <- which.max(m)
    states <- prod(m[-pivot] + 1)
    if (!is.finite(states) || states > control$max_states * (1 + 1e-12)) {
      stop("explicit groups require too many states; aggregate more indices or increase max_states")
    }

    upper_values <- lapply(groups, function(g) c(0, cumsum(weights[g])))
    lower_values <- lapply(groups, function(g) c(0, cumsum(rev(weights[g]))))
    errors <- vapply(groups, gap, numeric(1))
    deterministic_error <- sum(errors)
    exact_grouping <- all(vapply(groups, function(g) {
      all(weights[g] == weights[g[1L]])
    }, logical(1)))

    binomial_table <- function(n) {
      # Preserve exact dyadic probabilities for small fair-binomial groups.
      if (prob == 0.5 && n <= 52L) {
        row <- 1
        for (i in seq_len(n)) row <- c(row, 0) + c(0, row)
        mass <- row / 2^n
        cdf <- c(0, cumsum(mass))
        survival <- c(1, rev(cumsum(rev(mass[-1L]))), 0)
        return(list(mass=mass, log_mass=log(mass), cdf=cdf,
          survival=survival, log_cdf=log(cdf), log_survival=log(survival)))
      }
      k <- 0:n
      list(mass=stats::dbinom(k, n, prob),
        log_mass=stats::dbinom(k, n, prob, log=TRUE),
        cdf=stats::pbinom(c(-1, k), n, prob),
        survival=stats::pbinom(c(-1, k), n, prob, lower.tail=FALSE),
        log_cdf=stats::pbinom(c(-1, k), n, prob, log.p=TRUE),
        log_survival=stats::pbinom(c(-1, k), n, prob,
          lower.tail=FALSE, log.p=TRUE))
    }

    remainder_upper <- remainder_lower <- 0
    remainder_mass <- 1
    remainder_log_mass <- 0
    for (j in setdiff(seq_along(groups), pivot)) {
      tab <- binomial_table(m[j])
      remainder_upper <- as.vector(outer(remainder_upper, upper_values[[j]], "+"))
      remainder_lower <- as.vector(outer(remainder_lower, lower_values[[j]], "+"))
      remainder_mass <- as.vector(outer(remainder_mass, tab$mass, "*"))
      remainder_log_mass <- as.vector(outer(remainder_log_mass, tab$log_mass, "+"))
    }
    tab <- binomial_table(m[pivot])

    # Number of pivot support values at most x-r. Correct cancellation at the
    # comparison boundary using the actually represented sum r+value.
    support_index <- function(x, remainder, values) {
      index <- findInterval(x - remainder, values)
      repeat {
        good <- index > 0L
        bad <- which(good & (remainder + values[pmax(1L, index)] > x))
        if (!length(bad)) break
        index[bad] <- index[bad] - 1L
      }
      repeat {
        good <- index < length(values)
        bad <- which(good &
          (remainder + values[pmin(length(values), index+1L)] <= x))
        if (!length(bad)) break
        index[bad] <- index[bad] + 1L
      }
      index
    }
    log_sum <- function(x) {
      top <- max(x)
      if (!is.finite(top)) return(top)
      min(0, top + log(sum(exp(x-top))))
    }

    quantile_intervals <- function(remainder, values) {
      support_max <- max(remainder) + tail(values, 1)
      result <- matrix(NA_real_, nrow=length(p), ncol=2L,
        dimnames=list(NULL, c("lower", "upper")))
      resolved <- logical(length(p))
      target_ok <- function(x, i) {
        index <- support_index(x, remainder, values) + 1L
        if (use_lower[i]) {
          if (!is.null(numeric_lower) && log_lower[i] > -650) {
            return(sum(remainder_mass * tab$cdf[index]) >= numeric_lower[i])
          }
          log_sum(remainder_log_mass + tab$log_cdf[index]) >= log_lower[i]
        } else {
          if (!is.null(numeric_upper) && log_upper[i] > -650) {
            return(sum(remainder_mass * tab$survival[index]) <= numeric_upper[i])
          }
          log_sum(remainder_log_mass + tab$log_survival[index]) <= log_upper[i]
        }
      }
      neighbor <- function(x, above) {
        index <- support_index(x, remainder, values)
        if (above) {
          use <- index < length(values)
          if (!any(use)) return(Inf)
          min(remainder[use] + values[index[use]+1L])
        } else {
          use <- index > 0L
          if (!any(use)) return(-Inf)
          max(remainder[use] + values[index[use]])
        }
      }
      for (i in seq_along(p)) {
        if (log_lower[i] == -Inf || target_ok(0, i)) {
          result[i, ] <- 0
          resolved[i] <- TRUE
          next
        }
        if (log_upper[i] == -Inf) {
          result[i, ] <- support_max
          resolved[i] <- TRUE
          next
        }
        lo <- 0
        hi <- support_max
        for (iteration in seq_len(control$max_iterations)) {
          mid <- lo + (hi-lo) / 2
          if (mid <= lo || mid >= hi) break
          if (target_ok(mid, i)) hi <- mid else lo <- mid
        }
        # The successor of lo cannot exceed the true quantile; the predecessor
        # of hi has the same CDF as hi. Usually these identify one support point.
        lower <- neighbor(lo, TRUE)
        upper <- neighbor(hi, FALSE)
        if (is.finite(lower) && target_ok(lower, i)) {
          result[i, ] <- lower
          resolved[i] <- TRUE
        } else if (is.finite(upper) && upper >= lower && target_ok(upper, i)) {
          result[i, ] <- c(lower, upper)
          resolved[i] <- lower == upper
        } else {
          result[i, ] <- c(lo, hi)
        }
      }
      list(interval=result, resolved=resolved)
    }

    upper_search <- quantile_intervals(remainder_upper, upper_values[[pivot]])
    lower_search <- if (exact_grouping) upper_search else
      quantile_intervals(remainder_lower, lower_values[[pivot]])
    answer <- unname(pmin(total, upper_search$interval[, "upper"]))
    lower <- unname(pmin(total, pmax(0, lower_search$interval[, "lower"])))
    answer[log_lower == -Inf] <- lower[log_lower == -Inf] <- 0
    answer[log_upper == -Inf] <- lower[log_upper == -Inf] <- total
    if (!details) return(answer)
    search_error <- unname((upper_search$interval[, "upper"] -
      upper_search$interval[, "lower"]) +
      (lower_search$interval[, "upper"] - lower_search$interval[, "lower"]))
    group_summary <- data.frame(group=seq_along(groups), count=m,
      min_weight=vapply(groups, function(g) min(weights[g]), numeric(1)),
      max_weight=vapply(groups, function(g) max(weights[g]), numeric(1)),
      deterministic_error=errors)
    list(quantile=answer, lower=lower,
      deterministic_error=deterministic_error, groups=groups,
      group_summary=group_summary, pivot=pivot,
      enumerated_states=length(remainder_mass),
      total_count_states=prod(m+1), exact_grouping=exact_grouping,
      summary=data.frame(p=p, lower=lower, quantile=answer,
        bracket_width=answer-lower, search_error=search_error,
        resolved=upper_search$resolved & lower_search$resolved))
  }
  conditioned_quantiles <- function(
      p, weights, prob=0.5, groups=NULL, lower.tail=TRUE, log.p=FALSE,
      details=FALSE, lambda=NULL, max_states=100000,
      max_mgf_work=20000000, max_iterations=80) {
    flag <- function(x) is.logical(x) && length(x)==1L && !is.na(x)
    if (!flag(lower.tail) || !flag(log.p) || !flag(details))
      stop("lower.tail, log.p, and details must be single logical values.")
    if (!is.numeric(p) || anyNA(p) || any(p > if (log.p) 0 else 1) ||
        any(p < if (log.p) -Inf else 0) || any(is.nan(p)))
      stop("p must contain valid quantile probabilities.")
    if (!is.numeric(weights) || anyNA(weights) ||
        any(!is.finite(weights)) || any(weights < 0))
      stop("weights must contain finite nonnegative numbers.")
    if (!is.numeric(prob) || length(prob)!=1L || !is.finite(prob) ||
        prob < 0 || prob > 1)
      stop("prob must be one common Bernoulli probability in [0, 1].")
    for (z in list(max_states, max_mgf_work, max_iterations)) {
      if (!is.numeric(z) || length(z)!=1L || !is.finite(z) ||
          z < 1 || z != floor(z)) stop("Work limits must be positive integers.")
    }
    if (!is.null(lambda) && (!is.numeric(lambda) || anyNA(lambda) ||
        any(!is.finite(lambda)) || any(lambda < 0)))
      stop("lambda must be NULL or a finite nonnegative numerical vector.")
    n <- length(weights)
    if (!is.null(groups)) {
      if (!is.list(groups)) stop("groups must be a list of index vectors.")
      for (g in groups) {
        if (!is.numeric(g) || anyNA(g) || any(!is.finite(g)) ||
            any(g != floor(g)) || any(g < 1 | g > n))
          stop("Every group must contain valid integer indices into weights.")
      }
      used <- unlist(groups, use.names=FALSE)
      if (anyDuplicated(used)) stop("The groups must be disjoint.")
      groups <- c(groups, lapply(setdiff(seq_len(n), used), function(i) i))
    } else groups <- list(seq_len(n))
    groups <- lapply(groups, function(g) g[weights[g] > 0])
    groups <- Filter(length, groups)
    total <- sum(weights)
    if (!is.finite(total)) stop("The sum of weights must be finite.")
    log_complement <- function(x) {
      ans <- x
      low <- x < -log(2)
      ans[low] <- log1p(-exp(x[low]))
      ans[!low] <- log(-expm1(x[!low]))
      ans
    }
    lp <- if (log.p) p else log(p)
    log_alpha <- if (lower.tail) {
      if (log.p) log_complement(p) else log1p(-p)
    } else lp
    lower_endpoint <- if (lower.tail) lp == -Inf else lp == 0
    upper_endpoint <- if (lower.tail) lp == 0 else lp == -Inf
    probability <- if (lower.tail) exp(lp) else -expm1(lp)
    numeric_lower <- if (!log.p) { if (lower.tail) p else 1-p } else NULL
    numeric_upper <- if (!log.p) { if (lower.tail) 1-p else p } else NULL
    if (!length(groups) || prob %in% c(0, 1)) {
      value <- if (prob==1) total else 0
      out <- rep(value, length(p))
      if (!details) return(out)
      return(list(quantile=out, summary=data.frame(
        probability=probability, quantile=out, lower=out, prefix_upper=out,
        log_tail_bound=rep(-Inf, length(p)), exact=rep(TRUE,length(p))),
        groups=groups, count_states=1, deterministic_error=0,
        lambda=lambda, method="constant"))
    }
    sizes <- lengths(groups)
    log_states <- sum(log(sizes + 1))
    if (log_states > log(max_states) + 1e-12)
      stop("The joint count states exceed max_states; use fewer aggregated groups or retain fewer individual weights.")
    states <- prod(sizes + 1)
    if (!is.finite(states) || states > max_states)
      stop("The joint count states exceed max_states.")
    scale <- max(weights)
    w <- weights / scale
    lambda <- unique(c(0, lambda))
    use_mgf <- length(lambda) > 1L
    scaled_lambda <- lambda * scale
    if (any(!is.finite(scaled_lambda)))
      stop("The products of lambda and the largest weight must be finite.")
    if (use_mgf && (length(lambda)*sum(sizes*(sizes+1)/2) > max_mgf_work ||
                    length(lambda)*states > max_mgf_work))
      stop("The conditional MGF calculation exceeds max_mgf_work; use a smaller lambda grid or smaller groups.")
    log_add <- function(a, b) {
      z <- pmax(a,b)
      ans <- z
      ok <- is.finite(z)
      ans[ok] <- z[ok]+log1p(exp(-abs(a[ok]-b[ok])))
      ans
    }
    log_sum <- function(x) {
      z <- max(x)
      if (!is.finite(z)) return(z)
      z + log(sum(exp(x-z)))
    }
    # Conditional centered MGFs, one row for each count k=0,...,m.
    group_mgf <- function(a) {
      m <- length(a)
      z <- a-mean(a)
      ans <- matrix(0, nrow=m+1L, ncol=length(lambda))
      if (all(a==a[1L])) return(ans)
      for (ell in seq_along(lambda)) {
        if (lambda[ell]==0) next
        f <- 0
        for (j in seq_len(m)) {
          k <- 0:j
          no <- c(f, -Inf)+log((j-k)/j)
          yes <- c(-Inf, f)+scaled_lambda[ell]*z[j]+log(k/j)
          f <- log_add(no, yes)
        }
        ans[,ell] <- pmax(0,f)
      }
      ans[c(1L,m+1L),] <- 0
      ans
    }
    low <- high <- mu <- variance <- proxy <- log_mass <- 0
    mass <- 1
    log_mgf <- if (use_mgf) matrix(0, 1L, length(lambda)) else NULL
    deterministic_error <- 0
    for (g in groups) {
      a <- sort(w[g])
      m <- length(a)
      k <- 0:m
      a_low <- c(0,cumsum(a))
      a_high <- c(0,cumsum(rev(a)))
      a_mu <- k*mean(a)
      a_low[m+1L] <- a_high[m+1L] <- a_mu[m+1L] <- sum(a)
      a_var <- if (m > 1L) k*(m-k)/(m*(m-1))*sum((a-mean(a))^2) else c(0,0)
      # Bardenet--Maillard (2015), Propositions 2.2 and 2.3, combined.
      a_proxy <- pmin(k*(m-k+1), (k+1)*(m-k))/m*(max(a)-min(a))^2
      a_log <- stats::dbinom(k,m,prob,log=TRUE)
      if (prob==0.5 && m<=52L) {
        row <- 1
        for (j in seq_len(m)) row <- c(row,0)+c(0,row)
        a_mass <- row/2^m
        a_log <- log(a_mass)
      } else a_mass <- stats::dbinom(k,m,prob)
      previous <- length(mu)
      expand <- function(old, new) rep(old, each=m+1L)+rep(new, times=previous)
      low <- expand(low,a_low)
      high <- expand(high,a_high)
      mu <- expand(mu,a_mu)
      variance <- expand(variance,a_var)
      proxy <- expand(proxy,a_proxy)
      log_mass <- expand(log_mass,a_log)
      mass <- rep(mass,each=m+1L)*rep(a_mass,times=previous)
      if (use_mgf) {
        part <- group_mgf(a)
        log_mgf <- log_mgf[rep(seq_len(previous),each=m+1L),,drop=FALSE]+
          part[rep(seq_len(m+1L),times=previous),,drop=FALSE]
      }
      deterministic_error <- deterministic_error+max(a_high-a_low)
    }
    # CDF inversion for the two discrete count envelopes. We use sorted support
    # and a log-tail test, rather than subtracting a nearly unit CDF from one.
    envelope_quantile <- function(value, target, lower_bound=FALSE) {
      if (length(groups)==1L) {
        count <- stats::qbinom(p,sizes[1L],prob,lower.tail=lower.tail,log.p=log.p)
        return(value[count+1L])
      }
      ord <- order(value)
      s <- value[ord]
      log_weight <- log_mass[ord]
      weight <- mass[ord]
      cdf <- cumsum(weight)
      survival <- c(rev(cumsum(rev(weight[-1L]))),0)
      survivor <- rep(-Inf, length(s))
      running <- -Inf
      for (j in rev(seq_along(s))) {
        survivor[j] <- running
        running <- log_add(running,log_weight[j])
      }
      # A tied support value includes all its probability mass.
      last <- !duplicated(s,fromLast=TRUE)
      s <- s[last]
      survivor <- survivor[last]
      cdf <- cdf[last]
      survival <- survival[last]
      dyadic <- prob==0.5 && sum(sizes)<=52L
      guard <- if (lower_bound && !dyadic)
        64*.Machine$double.eps*(sum(sizes)+length(weight)) else 0
      vapply(seq_along(target),function(i) {
        if (target[i] >= -log(2)) {
          wanted <- if (is.null(numeric_lower))
            exp(log_complement(target[i])) else numeric_lower[i]
          wanted <- wanted/(1+guard)
          ok <- cdf >= wanted
        } else if (target[i] > -650) {
          wanted <- if (is.null(numeric_upper)) exp(target[i]) else numeric_upper[i]
          ok <- survival <= wanted*(1+guard)
        } else ok <- survivor <= target[i]+log1p(guard)
        if (any(ok)) s[which(ok)[1L]] else tail(s,1L)
      },numeric(1))
    }
    prefix <- envelope_quantile(high,log_alpha)
    lower <- envelope_quantile(low,log_alpha,lower_bound=TRUE)
    prefix[lower_endpoint] <- 0
    lower[lower_endpoint] <- 0
    prefix[upper_endpoint] <- sum(w)
    lower[upper_endpoint] <- sum(w)
    conditional_tail <- function(x) {
      log_bound <- rep(0,length(mu))
      log_bound[x >= high] <- -Inf
      constant <- variance==0
      log_bound[constant] <- ifelse(x < mu[constant],0,-Inf)
      active <- which(!constant & x < high & x > low)
      if (length(active)) {
        delta <- x-mu[active]
        spread <- high[active]-low[active]
        mean_above_low <- mu[active]-low[active]
        slack <- pmax(0,mean_above_low*(high[active]-mu[active])-variance[active])
        polynomial <- (mean_above_low+slack/(x-low[active]))/spread
        b <- pmin(0,log(polynomial))
        above <- which(delta > 0)
        if (length(above)) {
          lv <- log(variance[active[above]])
          cantelli <- lv-log_add(lv,2*log(delta[above]))
          b[above] <- pmin(b[above],cantelli)
          hp <- proxy[active[above]]
          use <- which(hp > 0)
          if (length(use)) b[above[use]] <- pmin(b[above[use]],
            -2*delta[above[use]]^2/hp[use])
        }
        log_bound[active] <- b
      }
      if (use_mgf) {
        for (ell in seq_along(lambda)) {
          b <- log_mgf[,ell]-scaled_lambda[ell]*(x-mu)
          log_bound <- pmin(log_bound,b)
        }
      }
      min(0,log_sum(log_mass+log_bound))
    }
    answer <- prefix
    log_tail <- vapply(prefix,conditional_tail,numeric(1))
    exact <- rep(all(variance==0),length(p))
    for (i in seq_along(p)) {
      if (lower_endpoint[i] || upper_endpoint[i] || exact[i]) next
      left <- 0
      right <- prefix[i]
      if (conditional_tail(left) <= log_alpha[i]) {
        answer[i] <- left
        next
      }
      # Keep an admissible upper bracket throughout. Stop at representable
      # resolution; the returned value need not coincide with a support point.
      for (iteration in seq_len(max_iterations)) {
        middle <- left+(right-left)/2
        if (middle==left || middle==right) break
        if (conditional_tail(middle) <= log_alpha[i]) right <- middle else left <- middle
      }
      answer[i] <- right
    }
    log_tail <- vapply(answer,conditional_tail,numeric(1))
    margin <- 128*.Machine$double.eps*max(1,n)*sum(w)
    answer <- pmin(sum(w),answer+ifelse(exact | lower_endpoint | upper_endpoint,0,margin))
    log_tail <- vapply(answer,conditional_tail,numeric(1))
    saturated <- answer >= sum(w)
    answer <- pmin(total,answer*scale)
    answer[saturated | upper_endpoint] <- total
    answer[lower_endpoint] <- 0
    lower <- pmax(0,lower-margin)*scale
    lower[lower_endpoint] <- 0
    lower[upper_endpoint] <- total
    if (!details) return(answer)
    list(quantile=answer, summary=data.frame(
      probability=probability, quantile=answer, lower=lower,
      prefix_upper=prefix*scale, log_tail_bound=log_tail, exact=exact),
      groups=groups, count_states=length(mu),
      deterministic_error=deterministic_error*scale,
      lambda=if (use_mgf) lambda else NULL,
      method=if (use_mgf) "conditional_moments_and_mgf" else "conditional_moments")
  }

  for (flag in list(lower.tail=lower.tail, log.p=log.p, details=details)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("lower.tail, log.p, and details must each be TRUE or FALSE.")
  }
  if (!is.numeric(p) || !is.null(dim(p)) || anyNA(p) ||
      (if (log.p) any(p > 0) else any(!is.finite(p) | p < 0 | p > 1)))
    stop("p must be a probability vector, or its logarithm when log.p=TRUE.")
  if (!is.numeric(weights) || !is.null(dim(weights)) ||
      any(!is.finite(weights) | weights < 0) || !is.finite(sum(weights)))
    stop("weights must be finite nonnegative numbers with a finite sum.")
  if (!is.numeric(prob) || !is.null(dim(prob)) ||
      !length(prob) %in% c(1L, length(weights)) ||
      any(!is.finite(prob) | prob < 0 | prob > 1))
    stop("prob must be a success probability, or one probability per weight.")
  if (!is.numeric(max_length) || length(max_length) != 1L ||
      !is.finite(max_length) || max_length < 0 || max_length > 52 ||
      max_length != floor(max_length))
    stop("max_length must be an integer between 0 and 52; the default is 25.")

  defaults <- list(max_states=2^18, max_dp_work=2e7,
    max_mitm_work=2e7, max_iterations=80, chernoff_iterations=60,
    max_p2_B=10000, max_p2_work=200000,
    max_group_states=65536, max_count_states=100000,
    max_mgf_work=2e7, mgf_points=37)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!names(control) %in% names(defaults)) ||
       anyDuplicated(names(control)))))
    stop("control must be a named list of documented resource limits.")
  defaults[names(control)] <- control
  control <- defaults
  optional <- c("max_p2_B", "max_p2_work", "max_group_states",
                "max_count_states", "max_mgf_work")
  for (name in names(control)) {
    value <- control[[name]]
    minimum <- if (name %in% optional) 0 else if (name == "mgf_points") 2 else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < minimum || value != floor(value))
      stop(paste0("Invalid control$", name, "."))
  }
  if (max(control$max_states, control$max_group_states,
          control$max_count_states) > 2^30)
    stop("State limits must be at most 2^30.")

  lower_prob <- upper_prob <- NULL
  if (log.p) {
    log_lower <- p
    log_upper <- log_complement(p)
  } else {
    lower_prob <- p
    upper_prob <- 1-p
    log_lower <- log(p)
    log_upper <- log1p(-p)
  }
  if (!lower.tail) {
    temporary <- log_lower
    log_lower <- log_upper
    log_upper <- temporary
    temporary <- lower_prob
    lower_prob <- upper_prob
    upper_prob <- temporary
  }
  n_request <- length(p)
  answer <- data.frame(probability=exp(log_lower),
    quantile=rep(NA_real_, n_request), exact=rep(FALSE, n_request),
    method=rep(NA_character_, n_request),
    tail_bound=probability_upper(log_upper), log_tail_bound=log_upper,
    grid_error=rep(NA_real_, n_request))
  bounds <- NULL
  diagnostics <- list(input_length=length(weights), active_length=0L,
    max_length=max_length, control=control, exact_attempts=list(),
    partitions=list())
  finish <- function() {
    q <- answer$quantile
    names(q) <- names(p)
    if (!details) return(q)
    list(quantile=q, summary=answer, bounds=bounds, diagnostics=diagnostics)
  }
  if (!n_request) return(finish())

  prob <- rep_len(prob, length(weights))
  offset <- sum(weights[prob == 1])
  active <- weights > 0 & prob > 0 & prob < 1
  w <- weights[active]
  pr <- prob[active]
  n <- length(w)
  diagnostics$active_length <- n
  if (!n) {
    answer$quantile <- offset
    answer$exact <- TRUE
    answer$method <- "constant"
    answer$tail_bound <- 0
    answer$log_tail_bound <- -Inf
    answer$grid_error <- 0
    return(finish())
  }
  order_w <- order(w)
  w <- w[order_w]
  pr <- pr[order_w]
  total <- sum(w)
  upper_endpoint <- offset+total
  round_margin <- 128*.Machine$double.eps*upper_endpoint+8*2^-1074
  upper_cap <- upper_endpoint+512*.Machine$double.eps*upper_endpoint+16*2^-1074
  if (!is.finite(upper_cap)) upper_cap <- Inf

  # Compare the smaller tail. In particular, never infer that q=0 merely
  # because 1-alpha and the zero-atom probability both rounded to one.
  log_zero <- sum(log1p(-pr))
  log_top <- sum(log(pr))
  fair <- all(pr == 0.5)
  zero_mass <- if (fair) 2^-n else exp(log_zero)
  top_mass <- if (fair) 2^-n else exp(log_top)
  endpoint_ok <- function(cdf, survival, log_cdf, log_survival) {
    vapply(seq_along(p), function(i) {
      if (log_lower[i] <= -log(2)) {
        if (!is.null(lower_prob) && log_lower[i] >= -650)
          return(cdf >= lower_prob[i])
        return(log_cdf >= log_lower[i])
      }
      if (!is.null(upper_prob) && log_upper[i] >= -650)
        return(survival <= upper_prob[i])
      log_survival <= log_upper[i]
    }, logical(1))
  }
  zero_survival <- if (fair && n <= 52L) 1-zero_mass else -expm1(log_zero)
  below_top <- if (fair && n <= 52L) 1-top_mass else -expm1(log_top)
  at_zero <- endpoint_ok(zero_mass, zero_survival,
    log_zero, log_complement(log_zero))
  at_top <- !endpoint_ok(below_top, top_mass,
    log_complement(log_top), log_top)
  at_zero <- at_zero | (log_lower == -Inf)
  at_top <- at_top | (log_upper == -Inf)
  at_zero[log_upper == -Inf] <- FALSE
  at_top[at_zero] <- FALSE
  answer$quantile[at_zero] <- offset
  answer$exact[at_zero] <- TRUE
  answer$method[at_zero] <- "lower_endpoint"
  answer$log_tail_bound[at_zero] <- log_complement(log_zero)
  answer$tail_bound[at_zero] <- probability_upper(log_complement(log_zero))
  answer$grid_error[at_zero] <- 0
  answer$quantile[at_top] <- upper_endpoint
  answer$exact[at_top] <- TRUE
  answer$method[at_top] <- "upper_endpoint"
  answer$log_tail_bound[at_top] <- -Inf
  answer$tail_bound[at_top] <- 0
  answer$grid_error[at_top] <- 0
  pending <- which(is.na(answer$quantile))
  short <- n <= max_length
  exact_control <- control
  exact_control$mitm_B <- max_length
  if (short && n <= 25L) {
    exact_control$max_states <- max(control$max_states, 2^ceiling(n/2))
    exact_control$max_mitm_work <- Inf
    exact_control$max_iterations <- max(2200, control$max_iterations)
  }
  accept_exact <- function(value, route) {
    available <- !is.null(value) && any(is.finite(value$quantile))
    diagnostics$exact_attempts[[route]] <<- list(available=available,
      reason=if (available) "Finite-distribution calculation succeeded." else
        "Structure, work limit, or a numerical CDF jump prevented this route.")
    if (!available) return(invisible(NULL))
    done <- which(is.finite(value$quantile))
    indices <- pending[done]
    answer$quantile[indices] <<- offset+value$quantile[done]
    answer$exact[indices] <<- TRUE
    answer$method[indices] <<- value$method[done]
    answer$grid_error[indices] <<- 0
    pending <<- which(is.na(answer$quantile))
    invisible(NULL)
  }
  exact_arguments <- function() list(
    log_lower=log_lower[pending], log_upper=log_upper[pending],
    control=exact_control,
    lower_prob=if (is.null(lower_prob)) NULL else lower_prob[pending],
    upper_prob=if (is.null(upper_prob)) NULL else upper_prob[pending])

  # Structured exact computations have priority at every input length.
  if (length(pending) && all(w == w[1L]) && all(pr == pr[1L])) {
    args <- c(list(n=n, step=w[1L], success=pr[1L], kind="binomial"),
              exact_arguments())
    accept_exact(do.call(exact_standard_quantiles, args), "binomial")
  }
  if (length(pending) && control$max_group_states > 0) {
    args <- c(list(w=w, prob=pr), exact_arguments())
    accept_exact(do.call(structured_count_quantiles, args), "binomial_blocks")
  }
  lattice <- if (length(pending)) exact_lattice(w) else NULL
  if (length(pending) && fair && n <= 52L && !is.null(lattice)) {
    args <- c(list(w_integer=lattice$weights, prob=pr), exact_arguments())
    value <- do.call(dp_quantiles, args)
    if (!is.null(value)) value$quantile <- lattice$step*value$quantile
    accept_exact(value, "lattice_dp")
  }
  if (length(pending) && fair && !is.null(lattice) &&
      identical(as.numeric(sort(lattice$weights)), as.numeric(seq_len(n)))) {
    states <- floor(n*(n+1)/4)+1
    if (states <= control$max_states && n*states <= control$max_dp_work) {
      args <- c(list(n=n, step=lattice$step, success=0.5, kind="wilcoxon"),
                exact_arguments())
      accept_exact(do.call(exact_standard_quantiles, args), "wilcoxon")
    }
  }
  if (length(pending) && !is.null(lattice)) {
    args <- c(list(w_integer=lattice$weights, prob=pr), exact_arguments())
    value <- do.call(dp_quantiles, args)
    if (!is.null(value)) value$quantile <- lattice$step*value$quantile
    accept_exact(value, "lattice_dp")
  }
  if (length(pending) && short) {
    args <- c(list(w=w, prob=pr), exact_arguments())
    accept_exact(do.call(exact_mitm_quantiles, args), "meet_in_the_middle")
    if (length(pending))
      stop("The required short-vector exact calculation did not resolve every ",
        "CDF jump within the exact work limits or double precision. ",
        "No conservative value has been substituted. Increase the exact work ",
        "limits for max_length > 25, or use higher precision for unresolved jumps.")
  }
  if (!length(pending)) return(finish())

  # Every remaining calculation supplies an upper bound. Normalization by a
  # binary unit avoids moment overflow. Tiny normalized terms are replaced by
  # deterministic maxima, preserving the direction of the bound.
  unit <- 2^max(-1074, min(1023, floor(log2(max(w)))))
  normalized <- w/unit
  tiny <- normalized < .Machine$double.xmin
  omitted <- sum(w[tiny])
  wn <- normalized[!tiny]
  pn <- pr[!tiny]
  nn <- length(wn)
  normalized_offset <- offset+omitted
  common <- all(pn == pn[1L])
  tab_rows <- list()
  add_rows <- function(method, quantile=NA_real_, log_tail=NA_real_,
                       exact=FALSE, reason="", partition="", states=NA_real_,
                       grid_error=NA_real_, original_units=FALSE) {
    q <- rep_len(quantile, length(pending))
    lt <- rep_len(log_tail, length(pending))
    good <- !is.na(q) & !is.na(lt) & lt <= log_upper[pending]
    if (!original_units)
      q <- pmin(upper_cap, normalized_offset+unit*q+round_margin)
    else q <- pmin(upper_cap, offset+q+round_margin)
    reason <- rep_len(reason, length(pending))
    reason[!good & reason == ""] <-
      "The numerical tail certificate was not resolved at the requested level."
    z <- data.frame(request=pending, probability=exp(log_lower[pending]),
      method=method, quantile=q, tail_bound=probability_upper(lt),
      log_tail_bound=lt, available=good, selected=FALSE,
      exact=rep(FALSE, length(pending)), reason=reason,
      partition=partition, states=rep_len(states, length(pending)),
      grid_error=rep_len(grid_error, length(pending)), stringsAsFactors=FALSE)
    tab_rows[[length(tab_rows)+1L]] <<- z
    invisible(NULL)
  }
  unique_targets <- unique(log_upper[pending])
  evaluated <- bound_candidates(wn, pn, unique_targets, control)
  evaluated <- evaluated[match(log_upper[pending], unique_targets)]
  old_names <- c("support", "endpoint", "chernoff", "hoeffding", "cantelli",
    "bernstein", "binomial_p2", "bentkus_gaussian", "bentkus_gaussian_k",
    "pinelis", "pinelis_k", "symmetric_chebyshev", "tomaszewski", "symmetry")
  rademacher_names <- c("bentkus_gaussian", "bentkus_gaussian_k", "pinelis",
    "pinelis_k", "symmetric_chebyshev", "tomaszewski", "symmetry")
  for (name in old_names) {
    q <- lt <- rep(NA_real_, length(pending))
    why <- rep("", length(pending))
    for (j in seq_along(pending)) {
      tab <- evaluated[[j]]$candidates
      rows <- which(tab$method == name)
      if (length(rows)) {
        r <- rows[which.min(ifelse(is.na(tab$quantile[rows]), Inf, tab$quantile[rows]))]
        q[j] <- tab$quantile[r]
        lt[j] <- tab$log_tail_bound[r]
      }
      if (is.na(q[j])) {
        why[j] <- if (name %in% rademacher_names && !all(pn == 0.5)) {
          "Requires fair active Bernoulli variables."
        } else if (name == "endpoint") {
          "The upper endpoint atom exceeds the requested tail probability."
        } else if (name == "tomaszewski") {
          "Requires an upper-tail probability at least 1/4."
        } else if (name == "symmetry") {
          "Requires an upper-tail probability at least 1/2."
        } else if (name %in% c("bentkus_gaussian_k", "pinelis_k")) {
          "The split coincides with, or does not improve, its unsplit comparison."
        } else if (name == "binomial_p2" &&
            (nn > control$max_p2_B || control$max_p2_work < nn+1)) {
          "Exceeds max_p2_B or max_p2_work."
        } else "A numerical certificate was not resolved within the work limits."
      }
    }
    add_rows(name, q, lt, reason=why)
  }

  # Upward rounding is evaluated alongside every concentration candidate.
  grid <- grid_quantiles(w, pr, log_lower[pending], log_upper[pending], control,
    lower_prob=if (is.null(lower_prob)) NULL else lower_prob[pending],
    upper_prob=if (is.null(upper_prob)) NULL else upper_prob[pending])
  if (is.null(grid)) add_rows("grid", reason="Exceeds max_states or max_dp_work.") else
    add_rows("grid", grid$quantile, log_upper[pending], grid$exact,
      grid_error=grid$error+round_margin, original_units=TRUE)

  # A fixed, disclosed collection of deterministic partitions is compared.
  # The greedy envelope partition adds further resolution under its budget.
  partitions <- list()
  partition_keys <- character()
  add_partition <- function(name, groups) {
    groups <- lapply(groups, function(g) sort(as.integer(g)))
    groups <- groups[lengths(groups) > 0L]
    groups <- groups[order(vapply(groups, min, numeric(1)))]
    key <- paste(vapply(groups, paste, collapse=",", FUN.VALUE=character(1)),
                 collapse=";")
    if (!key %in% partition_keys) {
      partitions[[name]] <<- groups
      partition_keys <<- c(partition_keys, key)
    }
  }
  grouped_cache <- list()
  if (common) {
    add_partition("one", list(seq_len(nn)))
    for (k in c(2L, 3L)) {
      if (nn >= k) {
        blocks <- split(seq_len(nn), pmin(k, floor((seq_len(nn)-1)*k/nn)+1L))
        add_partition(if (k == 2L) "two" else "three", blocks)
      }
    }
    for (r in c(1L, 2L, 4L)) {
      if (nn > r)
        add_partition(paste0("retain", r),
          c(list(seq_len(nn-r)), lapply(seq.int(nn-r+1L, nn), function(i) i)))
    }
    if (control$max_group_states > 0) {
      auto <- tryCatch(grouped_quantiles(
        p[pending], wn, prob=pn[1L], lower.tail=lower.tail, log.p=log.p,
        details=TRUE, control=list(max_states=control$max_group_states,
                                  max_iterations=control$max_iterations)),
        error=function(e) e)
      if (inherits(auto, "error")) {
        add_rows("grouped_auto", reason=conditionMessage(auto), partition="auto")
      } else {
        old_count <- length(partitions)
        add_partition("auto", auto$groups)
        if (length(partitions) > old_count) grouped_cache$auto <- auto
      }
    } else add_rows("grouped_auto", reason="Disabled by max_group_states=0.",
                    partition="auto")
  } else {
    add_rows("grouped", reason="Count envelopes require one common active success probability.")
    add_rows("conditional_moments", reason="Conditional count moments require one common active success probability.")
    add_rows("conditional_mgf", reason="Conditional count MGFs require one common active success probability.")
  }

  for (name in names(partitions)) {
    groups <- partitions[[name]]
    sizes <- lengths(groups)
    log_count <- sum(log1p(sizes))
    count_states <- if (log_count < log(.Machine$double.xmax)) exp(log_count) else Inf
    pivot_states <- prod(sizes[-which.max(sizes)]+1)
    count_states <- if (is.finite(count_states)) prod(sizes+1) else Inf
    meta <- list(sizes=sizes, pivot_states=pivot_states,
                 count_states=count_states, lambda=NULL)
    if (control$max_group_states == 0 || pivot_states > control$max_group_states) {
      add_rows(paste0("grouped_", name),
        reason="Exceeds max_group_states.", partition=name, states=pivot_states)
    } else {
      result <- grouped_cache[[name]]
      if (is.null(result)) result <- tryCatch(grouped_quantiles(
        p[pending], wn, prob=pn[1L], groups=groups,
        lower.tail=lower.tail, log.p=log.p, details=TRUE,
        control=list(max_states=control$max_group_states,
                     max_iterations=control$max_iterations)), error=function(e) e)
      if (inherits(result, "error"))
        add_rows(paste0("grouped_", name), reason=conditionMessage(result),
                 partition=name, states=pivot_states)
      else add_rows(paste0("grouped_", name), result$quantile,
        log_upper[pending], result$exact_grouping & result$summary$resolved,
        partition=name, states=pivot_states,
        grid_error=omitted+unit*(result$deterministic_error+result$summary$search_error)+round_margin)
    }
    if (control$max_count_states == 0 || count_states > control$max_count_states) {
      for (family in c("conditional_moments_", "conditional_mgf_"))
        add_rows(paste0(family, name), reason="Exceeds max_count_states.",
                 partition=name, states=count_states)
      diagnostics$partitions[[name]] <- meta
      next
    }
    # A slightly smaller tail target absorbs ordinary log-sum rounding.
    target <- log_upper[pending]-64*.Machine$double.eps*(1+abs(log_upper[pending]))
    run_conditioned <- function(lambda=NULL) tryCatch(conditioned_quantiles(
      target, wn, prob=pn[1L], groups=groups, lower.tail=FALSE,
      log.p=TRUE, details=TRUE, lambda=lambda,
      max_states=control$max_count_states,
      max_mgf_work=max(1, control$max_mgf_work),
      max_iterations=control$max_iterations), error=function(e) e)
    moments <- run_conditioned()
    if (inherits(moments, "error"))
      add_rows(paste0("conditional_moments_", name), reason=conditionMessage(moments),
               partition=name, states=count_states)
    else add_rows(paste0("conditional_moments_", name), moments$quantile,
      moments$summary$log_tail_bound, moments$summary$exact,
      partition=name, states=count_states)

    spread <- max(vapply(groups, function(g) diff(range(wn[g])), numeric(1)))
    mgf_work <- control$mgf_points*max(sum(sizes*(sizes+1)/2), count_states)
    if (control$max_mgf_work == 0 || mgf_work > control$max_mgf_work) {
      add_rows(paste0("conditional_mgf_", name),
        reason="Exceeds max_mgf_work.", partition=name, states=count_states)
    } else if (spread == 0) {
      if (inherits(moments, "error"))
        add_rows(paste0("conditional_mgf_", name), reason=conditionMessage(moments),
          partition=name, states=count_states)
      else add_rows(paste0("conditional_mgf_", name), moments$quantile,
        moments$summary$log_tail_bound, moments$summary$exact,
        reason="Counts determine each group exactly; its conditional MGF is constant.",
        partition=name, states=count_states)
    } else {
      lambda <- c(0, exp(seq(log(0.01), log(500),
                             length.out=control$mgf_points-1L))/spread)
      meta$lambda <- lambda/unit
      if (any(!is.finite(lambda)))
        add_rows(paste0("conditional_mgf_", name),
          reason="The automatic lambda grid exceeds floating-point range.",
          partition=name, states=count_states)
      else {
        mgf <- run_conditioned(lambda)
        if (inherits(mgf, "error"))
          add_rows(paste0("conditional_mgf_", name), reason=conditionMessage(mgf),
            partition=name, states=count_states)
        else add_rows(paste0("conditional_mgf_", name), mgf$quantile,
          mgf$summary$log_tail_bound, mgf$summary$exact,
          partition=name, states=count_states)
      }
    }
    diagnostics$partitions[[name]] <- meta
  }
  bounds <- do.call(rbind, tab_rows)
  rownames(bounds) <- NULL
  for (i in pending) {
    rows <- which(bounds$request == i & bounds$available)
    if (!length(rows)) stop("No valid quantile candidate was available.")
    winner <- rows[which.min(bounds$quantile[rows])]
    bounds$selected[winner] <- TRUE
    answer$quantile[i] <- bounds$quantile[winner]
    answer$exact[i] <- bounds$exact[winner]
    answer$method[i] <- bounds$method[winner]
    answer$tail_bound[i] <- bounds$tail_bound[winner]
    answer$log_tail_bound[i] <- bounds$log_tail_bound[winner]
    answer$grid_error[i] <- bounds$grid_error[winner]
  }
  finish()
}
