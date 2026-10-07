\name{qweighted_bernoulli}
\alias{qweighted_bernoulli}
\title{Ordinary and Randomized Quantiles of Weighted Bernoulli Sums}
\description{
Computes exact finite-distribution quantiles when the number of random
summands is small or an affordable exact representation is available.
Otherwise evaluates the applicable conservative quantile candidates within
the work limits and returns their minimum. All methods are accessed
through this one function. Optional randomization returns the adjacent
cutoffs and their mixing probability, including the original inclusive
two-sided construction for fair Bernoulli sums.
}
\usage{
qweighted_bernoulli(
  p, weights, prob = 0.5, max_length = 25,
  lower.tail = TRUE, log.p = FALSE,
  details = FALSE, control = list(),
  randomized = FALSE, two.sided = FALSE, draw = TRUE
)
}
\arguments{
  \item{p}{Numeric vector of quantile probabilities. Values must be in
  \code{[0, 1]}, or in \code{[-Inf, 0]} when \code{log.p = TRUE}.
  With \code{lower.tail = FALSE}, supplies upper-tail probabilities.
  Names are preserved on the returned quantile vector. Missing values
  are not allowed; an empty vector is allowed.}

  \item{weights}{Numeric vector of finite, nonnegative weights, with a
  finite sum. Zero weights and an empty vector are allowed.}

  \item{prob}{Bernoulli success probabilities: either one number recycled
  to every weight, or a numeric vector with the same length as
  \code{weights}. Values must be finite and in \code{[0, 1]}.
  These are the probabilities of the independent Bernoulli variables,
  distinct from the requested quantile levels in \code{p}.}

  \item{max_length}{Integer from 0 to 52, default 25. Generic exact
  subset-sum calculation is required when the number of active summands
  is at most this value. An active summand has a positive weight and a
  success probability strictly between zero and one. Affordable exact
  structural representations are checked at any length. See
  \sQuote{Automatic calculation} and \sQuote{Work controls}.}

  \item{lower.tail}{Logical scalar. If \code{TRUE}, \code{p} specifies
  lower-tail quantile levels. If \code{FALSE}, \code{p} specifies
  upper-tail levels, so an ordinary-scale input \eqn{\alpha} requests
  the \eqn{1-\alpha} quantile.}

  \item{log.p}{Logical scalar. If \code{TRUE}, \code{p} contains the
  natural logarithms of the probabilities specified by \code{lower.tail}.}

  \item{details}{Logical scalar. In ordinary mode, if \code{FALSE}, return
  the numeric quantile vector; if \code{TRUE}, also return a summary, the
  candidate table, and calculation diagnostics. Randomized mode always
  returns a list with the mixing distribution and summary; \code{TRUE}
  adds the candidate table and calculation diagnostics.}

  \item{control}{A named list overriding any of the 11 work controls
  described below. Names must be distinct and recognized. Omitted
  entries use their defaults. These controls limit computational effort;
  they do not select a separate public method.}

  \item{randomized}{Logical scalar. If \code{TRUE}, return a mixing
  distribution for a randomized cutoff and, by default, an independent
  draw. The default \code{FALSE} retains ordinary quantile computation.}

  \item{two.sided}{Logical scalar. If \code{TRUE}, compute the original
  inclusive two-sided cutoff for a fair Bernoulli sum. Requires
  \code{randomized = TRUE}, \code{prob = 0.5} on every positive weight,
  and a requested tail level strictly between zero and one.}

  \item{draw}{Logical scalar. In randomized mode, if \code{TRUE}, generate
  a cutoff with the returned mixing distribution. If \code{FALSE}, return
  the distribution with \code{kappa = NULL} and \code{quantile = NULL}
  without consuming random numbers. Ordinary calls and degenerate
  mixtures also consume no random numbers.}
}
\details{
Let
\deqn{G = \sum_{i=1}^B w_i \xi_i, \qquad
       \xi_i \sim \mathrm{Bernoulli}(\pi_i),}{
       G = sum(w[i] * xi[i]), with xi[i] ~ Bernoulli(prob[i]),}
where the Bernoulli variables are independent. For \eqn{0 < u < 1},
the quantile convention is
\deqn{q_G(u) = \inf\{x : P(G \le x) \ge u\}.}{
       q_G(u) = inf\{x : P(G <= x) >= u\}.}
In ordinary mode, levels zero and one return the lower and upper support
endpoints, respectively. Deterministic summands are retained as an offset.

An exact answer computes this finite-distribution quantile. A conservative
answer \eqn{\widehat q} is an upper bound for it. In particular, a request
with upper-tail level \eqn{\alpha} is based on the \emph{strict-tail}
condition
\deqn{P(G > \widehat q) \le \alpha.}{P(G > qhat) <= alpha.}
A conservative answer need not be a support point. This condition is
different from one with \eqn{G \ge \widehat q}; an inclusive rejection rule
must account for a possible atom at the threshold.

Pass small upper-tail levels directly with \code{lower.tail = FALSE},
preferably also \code{log.p = TRUE}. Forming \code{1 - alpha} first can
round the input to one and discard the intended tail level. Multiple
levels in one call share distribution and bound preparation.
}
\section{Randomized cutoffs}{
The tail level is \eqn{\alpha = 1-p} when \code{lower.tail = TRUE}, or
\eqn{\alpha = p} when \code{lower.tail = FALSE}, after interpreting
\code{log.p}. The following formulas describe the exact branch and
\eqn{0 < \alpha < 1}. The mixing probability \code{tau} always selects
the larger cutoff \code{c_plus}.

For one-sided randomization, \eqn{c_+} is the ordinary
\eqn{1-\alpha} quantile and \eqn{c_-} is its support predecessor,
with \eqn{-\infty} below the minimum support point. Define
\deqn{P_+ = P(G > c_+),\quad P_- = P(G > c_-),\quad
       \tau = \frac{P_- - \alpha}{P_- - P_+}.}{
       P_plus = P(G > c_plus), P_minus = P(G > c_minus),
       tau = (P_minus - alpha) / (P_minus - P_plus).}
An independent mixture choosing \eqn{c_+} with probability \eqn{\tau}
has \eqn{P(G > \kappa)=\alpha}. One-sided endpoint requests are allowed:
ordinary lower-tail level zero gives \eqn{-\infty}, and level one gives
the maximum support point.

For \code{two.sided = TRUE}, put \eqn{N=\sum_i w_i} and
\eqn{\mathcal S=\{-1\}\cup(\mathrm{supp}(G)\cap[0,N/2])}.
For \eqn{c\in\mathcal S} define
\deqn{P_h(c)=P(|G-N/2|\ge N/2-c),\qquad P_h(-1)=0.}{
       P_h(c) = P(abs(G - N/2) >= N/2 - c), with P_h(-1) = 0.}
Then \eqn{c_+=q_G(\alpha/2)} is the first support point where
\eqn{P_h(c)\ge\alpha}, and \eqn{c_-} is its predecessor in
\eqn{\mathcal S}. In this mode the probability fields are
\eqn{P_+=P_h(c_+)} and \eqn{P_-=P_h(c_-)} and
\deqn{\tau=\frac{\alpha-P_-}{P_+-P_-}.}{
       tau = (alpha - P_minus) / (P_plus - P_minus).}
The independent mixture satisfies
\eqn{P(|G-N/2|\ge N/2-\kappa)=\alpha}. Symmetry gives
\eqn{P_h(c)=2P(G\le c)} only for \eqn{c<N/2}; at the center,
\eqn{P_h(N/2)=1}. The algorithm counts a center atom once. For all-zero
weights it returns \eqn{c_+=0}, \eqn{c_-=-1}, and \eqn{\tau=\alpha}.

Exact engines return the boundary atom and actual predecessor from the
same law used to calculate the quantile. A predecessor need not be one
unit or one lattice step below the quantile. Exact randomization does
not require full enumeration of the Cartesian product in a
meet-in-the-middle calculation.

If the exact mixing distribution cannot be obtained for a longer
unstructured sum, the function returns a conservative degenerate mixture:
\code{c_plus = c_minus}, \code{tau = 1}, and \code{exact = FALSE}.
Its probability fields are upper bounds, not exact event probabilities,
and it does not assert an exact adjacent bracket or exact level attainment.
One-sided mode uses the minimum available upper quantile bound. Two-sided
mode first selects the minimum available strict upper quantile bound
\eqn{q} at level \eqn{\alpha/2}, then uses ascending and descending
greedy subset searches to find an actual support value \eqn{c<N-q},
or the sentinel \eqn{-1}. Symmetry implies
\eqn{P_h(c)\le 2P(G>q)\le\alpha}. The strict support comparison is
essential for the inclusive event. Selection of the deterministic bound
precedes any random draw.

The probability identities include a fresh independent uniform randomizer.
Do not reuse a data-dependent or tie-breaking uniform from the statistic.
With \code{draw = FALSE}, the mixing distribution is returned without
changing the RNG state. The mathematical note's randomization addendum
contains the proofs and further examples.
}
\section{Automatic calculation}{
The function first removes zero and deterministic terms, handles
endpoints, and checks affordable exact structure. The structural routes
include an equal-weight binomial distribution, repeated exact
\code{(weight, probability)} pairs represented by independent binomial
counts, a verified represented lattice, and fair rank weights. Structure
is recognized from the supplied floating-point values without approximate
rounding of weights to declare them equal or lattice-valued.

For unresolved requests with at most \code{max_length} active summands,
an exact meet-in-the-middle calculation is used. With the default
\code{max_length = 25}, a general length-25 problem has half-distributions
of 4096 and 8192 states. If a required exact calculation cannot resolve
a CDF jump within the applicable limits or double precision, the function
raises an error; it does not silently substitute a conservative answer.

Above \code{max_length}, unresolved requests proceed directly to the
conservative calculations, without generic subset-sum enumeration.
Candidates include support bounds, upward-grid quantiles, Chernoff,
Hoeffding, Cantelli, Bernstein, and binomial \eqn{P_2} comparisons.
Fair Bernoulli variables additionally permit Bentkus--Dzindzalieta and
Pinelis Gaussian comparisons, their Montgomery-Smith \eqn{K}-functional
refinements, and symmetry bounds. Grouped count envelopes and conditional
moment, support, Serfling, and finite-grid MGF refinements are also
compared when applicable.

The count-based conservative refinements currently require one common
success probability for the remaining random sum. The other applicable
candidates remain available for unequal probabilities. The count
calculations consider one full group, two or three blocks of ordered
weights, partitions retaining one, two, or four largest weights as
singletons, and an affordable greedy refinement; duplicate partitions
are removed. For largest positive within-group spread \eqn{R}, the
automatic conditional MGF grid consists of zero and
\code{mgf_points - 1} logarithmically spaced parameters from
\eqn{0.01/R} to \eqn{500/R}.

The returned bound is the minimum of the implemented, applicable
candidates completed within the work limits. It is not an optimization
over all possible partitions, all MGF parameters, or all probability
inequalities. Exact answers bypass unnecessary conservative work.
}
\section{Work controls}{
The recognized entries of \code{control} are:
\describe{
  \item{\code{max_states = 2^18}}{Maximum lattice/grid array length or
  meet-in-the-middle half-array length. The required exact route with at
  most 25 active terms raises the latter limit when necessary.}
  \item{\code{max_dp_work = 2e7}}{Maximum cumulative lattice/grid array
  updates for one PMF calculation; also determines an affordable upward
  grid.}
  \item{\code{max_mitm_work = 2e7}}{Maximum number of visited left-half
  entries in generic exact support queries. This cap is not imposed for
  a required exact calculation with at most 25 active terms.}
  \item{\code{max_iterations = 80}}{Iteration budget for scalar searches,
  bound inversions, and structured exact inversions. For a required exact
  calculation with at most 25 active terms, the exact branch uses at
  least 2200 support-search iterations.}
  \item{\code{chernoff_iterations = 60}}{MGF optimization evaluation
  budget per distinct tail level for the full Chernoff candidate.}
  \item{\code{max_p2_B = 10000}}{Maximum active count for the optional
  binomial \eqn{P_2} candidate. Zero disables this candidate.}
  \item{\code{max_p2_work = 200000}}{Preparation and scan budget for the
  binomial \eqn{P_2} calculation, shared across the requested levels.
  Zero disables this candidate.}
  \item{\code{max_group_states = 65536}}{Maximum enumerated non-pivot
  states for grouped count envelopes or exact repeated-class
  representations. The largest count is integrated through its binomial
  CDF. Zero disables these optional routes.}
  \item{\code{max_count_states = 100000}}{Maximum full joint count states
  for conditional moment and MGF calculations. Zero disables these
  conditional calculations.}
  \item{\code{max_mgf_work = 2e7}}{Conditional MGF preparation and
  joint-table work budget, applied separately to each candidate partition.
  Zero disables conditional MGF calculations while retaining applicable
  conditional moment calculations.}
  \item{\code{mgf_points = 37}}{Total number of automatic conditional
  MGF parameters, including zero; must be at least two. The default uses
  zero and 36 positive parameters.}
}
All controls must be finite scalar integers. The optional
\code{max_p2_B}, \code{max_p2_work}, \code{max_group_states},
\code{max_count_states}, and \code{max_mgf_work} controls may be zero.
Other controls must be positive, with \code{mgf_points >= 2}.
The three state caps \code{max_states}, \code{max_group_states}, and
\code{max_count_states} must each be at most \code{2^30}.

For a required exact calculation with at most 25 active terms, the
implementation permits at least \eqn{2^{\lceil n/2\rceil}} entries per
half-array, bypasses the query-work cap, and increases the support-search
budget as described above. Thus lowering conservative work budgets does
not change such a request into an approximation. Raising
\code{max_length} above 25 requests potentially larger exact calculations;
these are subject to the specified work limits and can raise an error if
unresolved. Setting \code{max_length = 0} disables only generic exact
subset-sum calculation, not structural exact routes.

The work limits are operation and state-count proxies, not elapsed-time
or memory-allocation guarantees. Automatic partitions or grids may
change when limits are changed, so a larger budget does not promise a
monotone decrease in the returned value.
}
\value{
For \code{randomized = FALSE} and \code{details = FALSE}, a numeric
vector of length \code{length(p)},
with the names of \code{p}.

For \code{randomized = FALSE} and \code{details = TRUE}, a list with
four components:
\describe{
  \item{\code{quantile}}{The same numeric quantile vector.}
  \item{\code{summary}}{A data frame with one row per requested level.
  Its columns are \code{probability} (the ordinary-scale lower-tail
  level), \code{quantile}, \code{exact}, \code{method},
  \code{tail_bound}, \code{log_tail_bound}, and \code{grid_error}.
  \code{method} identifies the selected internal engine, and
  \code{exact} identifies an exact finite-distribution calculation.
  Tail fields are upper bounds for the strict tail at the returned
  threshold, not measurements of its achieved probability; an exact
  interior answer usually reports the requested upper-tail level.
  \code{grid_error} is zero for exact answers, an absolute allowance
  for selected grid or grouped-envelope candidates, and \code{NA} for
  other selected bounds.}
  \item{\code{bounds}}{A data frame containing candidate rows for
  requests that entered the conservative branch, or \code{NULL} if all
  requests were resolved exactly. Columns are \code{request} (index
  into \code{p}), \code{probability}, \code{method}, \code{quantile},
  \code{tail_bound}, \code{log_tail_bound}, \code{available},
  \code{selected}, \code{exact}, \code{reason}, \code{partition},
  \code{states}, and \code{grid_error}. Unavailable candidates retain
  an explanation of the failed hypothesis, resource limit, or numerical
  calculation. For each request in this table, the selected quantile is
  the minimum over its available rows.}
  \item{\code{diagnostics}}{A list containing \code{input_length},
  \code{active_length}, \code{max_length}, the completed \code{control}
  list, \code{exact_attempts}, and \code{partitions}. Partition records
  contain group sizes, enumerated pivot-remainder and full count-state
  counts, and any MGF grid in the original weight units.}
}
The ordinary-scale \code{probability} field can round to one for a very
small upper tail; the logarithmic tail fields preserve the supplied
logarithmic level. Detailed fields describe the calculation and are not
an independent certification of floating-point error.

For \code{randomized = TRUE}, the return value is always a list:
\describe{
  \item{\code{c_plus}, \code{c_minus}}{Exact adjacent cutoffs or equal
  endpoints for a degenerate conservative fallback. Names are preserved.}
  \item{\code{tau}}{Probability of selecting \code{c_plus}.}
  \item{\code{log_tau}, \code{log_one_minus_tau}}{Natural logarithms of
  the two mixing probabilities, retained when ordinary-scale values
  round to zero or one.}
  \item{\code{P_plus}, \code{P_minus}}{The event probabilities at the
  cutoffs in the exact branch, or upper bounds in a conservative row.
  For one-sided mode the event is \eqn{G>c}; for two-sided mode it is
  \eqn{|G-N/2|\ge N/2-c}.}
  \item{\code{log_P_plus}, \code{log_P_minus}}{Natural logarithms of
  the corresponding probability fields.}
  \item{\code{kappa}, \code{quantile}}{The same generated cutoff vector.
  Both fields are \code{NULL} when \code{draw = FALSE}.}
  \item{\code{exact}}{One logical value per request, indicating that
  the exact-law mixing distribution was resolved.}
  \item{\code{alpha}, \code{log_alpha}, \code{two.sided}}{Requested
  tail levels, their logarithms, and the event convention.}
  \item{\code{summary}}{A data frame with columns \code{probability},
  \code{alpha}, \code{c_minus}, \code{c_plus}, \code{tau},
  \code{P_minus}, \code{P_plus}, \code{exact}, \code{degenerate},
  \code{probability_type}, and \code{method}. The probability type is
  \code{"exact"} or \code{"upper_bound"}.}
  \item{\code{bounds}, \code{diagnostics}}{Added when
  \code{details = TRUE}. In conservative two-sided mode, the candidate
  quantiles concern the underlying upper-tail level \eqn{\alpha/2},
  before reflection and the support search; they are not the final
  two-sided cutoffs.}
}
}
\note{
All proofs concern real arithmetic. Here \dQuote{exact} means an exact
finite-distribution algorithm evaluated in double precision, not symbolic
or certified interval arithmetic. Numerically indistinguishable subset
sums and probabilities at a CDF jump retain this qualification.
A randomized calculation raises an error if the adjacent support values
collapse to the same representable double after an offset is added.

Small tails are compared directly and logarithmic probabilities are
used when ordinary masses can underflow. Conservative calculations use
scaling, directional rounding safeguards, and feasible upper brackets.
Extremely small terms in a normalized conservative calculation may be
replaced by their deterministic maximum contribution, which is added
back to the bound. These safeguards do not constitute a complete
outward error analysis of R's special functions. The accompanying
mathematical note gives the derivations and numerical qualifications.

Draws use \code{\link[stats]{runif}} and the currently selected R random
number generator. Its finite resolution entails a small discretization
error relative to an ideal continuous uniform. Extremely small mixing
probabilities can lie below that resolution even when their logarithms
are retained. Use \code{draw = FALSE} to obtain the distribution for a
different sampling implementation if necessary; see \code{\link[base]{RNG}}.
The installed baseline note documents version 0.1.0. The current source
and \file{randomized_quantiles.pdf} addendum document version 0.2.0;
the baseline note's code appendix is not the current implementation.
}
\references{
Bentkus, V. (2004). On Hoeffding's inequalities.
\emph{The Annals of Probability}, \bold{32}(2), 1650--1673.
\doi{10.1214/009117904000000360}.

Bentkus, V. K. and Dzindzalieta, D. (2015). A tight Gaussian bound for
weighted sums of Rademacher random variables.
\emph{Bernoulli}, \bold{21}(2), 1231--1237.
\doi{10.3150/14-BEJ603}.

Pinelis, I. (2012). An asymptotically Gaussian bound on the Rademacher
tails. \emph{Electronic Journal of Probability}, \bold{17}, article 35.
\doi{10.1214/EJP.v17-2026}.

Montgomery-Smith, S. J. (1990). The distribution of Rademacher sums.
\emph{Proceedings of the American Mathematical Society},
\bold{109}(2), 517--522.
\doi{10.1090/S0002-9939-1990-1013975-0}.

Bardenet, R. and Maillard, O.-A. (2015). Concentration inequalities for
sampling without replacement. \emph{Bernoulli}, \bold{21}(3), 1361--1385.
\doi{10.3150/14-BEJ605}.

Keller, N. and Klein, O. (2022). Proof of Tomaszewski's conjecture on
randomly signed sums. \emph{Advances in Mathematics}, \bold{407},
article 108558. \doi{10.1016/j.aim.2022.108558}.
}
\seealso{
\code{\link[stats]{qbinom}}, \code{\link[stats]{qsignrank}},
\code{\link[stats]{runif}}, \code{\link[base]{RNG}},
\code{\link{weightedBernoulli-package}}.
}
\examples{
# Original fair two-sided mixing distribution, without an RNG draw.
cut <- qweighted_bernoulli(
  0.05, weights = 1:20, lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE, draw = FALSE
)
cut[c("c_minus", "c_plus", "tau", "P_minus", "P_plus")]
cut$tau * cut$P_plus + (1 - cut$tau) * cut$P_minus

# An atom at the center is counted once.
center <- qweighted_bernoulli(
  0.75, weights = c(1, 1), lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE, draw = FALSE
)
center$tau  # 0.5

# Reproducible draw from the exact mixing distribution.
set.seed(20261007)
qweighted_bernoulli(
  0.05, weights = 1:20, lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE
)$kappa

# One-sided randomization permits unequal Bernoulli probabilities.
one <- qweighted_bernoulli(
  0.05, weights = c(1, 2, 4), prob = c(0.2, 0.5, 0.8),
  lower.tail = FALSE, randomized = TRUE, draw = FALSE
)
one[c("c_minus", "c_plus", "tau")]

# Small vector with unequal probabilities: exact.
qweighted_bernoulli(
  c(0.5, 0.9, 0.95), weights = c(1, 2, 4),
  prob = c(0.2, 0.5, 0.8)
)

# Repeated weights give an exact calculation above length 25.
fit <- qweighted_bernoulli(
  0.05, weights = c(60, rep(1, 39)),
  lower.tail = FALSE, details = TRUE
)
fit$quantile
fit$summary[, c("quantile", "exact", "method")]

# A thousand equal weights reduce to one binomial distribution.
qweighted_bernoulli(0.95, weights = rep(1, 1000))

\donttest{
# Longer unstructured vector: compare all affordable candidates.
fit <- qweighted_bernoulli(
  c(0.05, 0.01), weights = sqrt(1:100),
  lower.tail = FALSE, details = TRUE
)
fit$summary
tab <- subset(fit$bounds, request == 1 & available)
tab <- tab[order(tab$quantile), ]
tab[, c("method", "quantile", "selected")]
stopifnot(fit$quantile[1] == min(tab$quantile))

# Supply very small upper tails without subtracting from one.
qweighted_bernoulli(
  log(1e-30), weights = sqrt(1:100),
  lower.tail = FALSE, log.p = TRUE
)
}
}
\keyword{distribution}
