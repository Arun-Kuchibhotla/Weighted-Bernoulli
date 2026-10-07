\name{weightedBernoulli-package}
\alias{weightedBernoulli-package}
\alias{weightedBernoulli}
\docType{package}
\title{Exact, Randomized, and Conservative Weighted Bernoulli Quantiles}
\description{
Quantiles for sums of independent Bernoulli random variables with finite,
nonnegative weights and scalar or individual success probabilities.
One exported function, \code{\link{qweighted_bernoulli}}, combines exact
distribution calculations, randomized cutoffs, and conservative
probability comparisons.
}
\details{
The default \code{max_length = 25} requests an exact finite-distribution
quantile for at most 25 active summands, after removing zero weights and
deterministic terms. Affordable exact structure is recognized at any
length, including equal weights, repeated weight--probability classes,
fair rank weights, and verified lattices. Longer unstructured problems
use the minimum of the applicable conservative candidates completed
within the documented work limits.

The conservative collection includes numerical grid and binomial-count
aggregation bounds, conditional refinements, and concentration and
Gaussian comparisons. There is no public method selector and no need
to call different functions to compare methods. Detailed output exposes
the selected calculation, candidate values, and reasons for unavailable
candidates.

Ordinary calls have the upper-tail contract \eqn{P(G > q) \le \alpha},
with a strict inequality in the event. Version 0.2.0 adds
\code{randomized = TRUE}, returning the adjacent exact cutoffs,
their mixing probability \code{tau}, the boundary event probabilities,
and an optional independent draw \code{kappa}. Set \code{draw = FALSE}
to obtain this distribution without consuming random numbers.

One-sided randomization calibrates \eqn{P(G>\kappa)=\alpha} in the exact
branch. With \code{two.sided = TRUE}, fair Bernoulli sums instead use the
original inclusive event \eqn{|G-N/2|\ge N/2-\kappa}, where
\eqn{N=\sum_i w_i}, with the sentinel \eqn{-1} below the support and
correct treatment of a center atom. In both cases \code{tau} selects
the larger cutoff. If exact computation is unaffordable, the function
returns a degenerate conservative cutoff marked \code{exact = FALSE};
its probability fields are upper bounds and do not assert exact nominal
attainment.

All computations use double precision;
\dQuote{exact} describes the finite-distribution algorithm and is not
a claim of symbolic or certified interval arithmetic.

The package uses base R and the standard \pkg{stats} and \pkg{utils} packages. The repository includes
installation and usage instructions, mathematical notes with proofs
and citations, reproducible comparisons, and validation scripts.
The version 0.2.0 randomization addendum can be located with
\code{system.file("doc", "randomized_quantiles.pdf", package = "weightedBernoulli")}.
The main note, \file{weighted_bernoulli_quantiles.pdf}, remains the
version 0.1.0 baseline; its printed code appendix is not the current
implementation. The package source and randomization addendum describe
the current interface. Draws use R's finite-resolution uniform generator;
the exact probability identities refer to an ideal independent uniform.
See \code{\link{qweighted_bernoulli}} for the complete interface,
work controls, numerical qualifications, and examples.
}
\seealso{\code{\link{qweighted_bernoulli}}}
\examples{
cut <- qweighted_bernoulli(
  0.05, weights = 1:20, lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE, draw = FALSE
)
cut[c("c_minus", "c_plus", "tau")]

qweighted_bernoulli(
  c(0.5, 0.9, 0.95), weights = c(1, 2, 4),
  prob = c(0.2, 0.5, 0.8)
)
}
\keyword{package}
\keyword{distribution}
