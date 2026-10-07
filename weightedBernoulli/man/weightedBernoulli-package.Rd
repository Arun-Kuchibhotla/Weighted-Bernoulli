\name{weightedBernoulli-package}
\alias{weightedBernoulli-package}
\alias{weightedBernoulli}
\docType{package}
\title{Exact and Conservative Quantiles of Weighted Bernoulli Sums}
\description{
Quantiles for sums of independent Bernoulli random variables with finite,
nonnegative weights and scalar or individual success probabilities.
One exported function, \code{\link{qweighted_bernoulli}}, combines exact
distribution calculations and conservative probability comparisons.
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

The upper-tail contract is \eqn{P(G > q) \le \alpha}, with a strict
inequality in the event. All computations use double precision;
\dQuote{exact} describes the finite-distribution algorithm and is not
a claim of symbolic or certified interval arithmetic.

The package uses base R and the standard \pkg{stats} and \pkg{utils} packages. The repository includes
installation and usage instructions, a mathematical note with proofs
and citations, reproducible comparisons, and independent validation
scripts. The installed note can be located with
\code{system.file("doc", "weighted_bernoulli_quantiles.pdf", package = "weightedBernoulli")}.
See \code{\link{qweighted_bernoulli}} for the complete interface,
work controls, numerical qualifications, and examples.
}
\seealso{\code{\link{qweighted_bernoulli}}}
\examples{
qweighted_bernoulli(
  c(0.5, 0.9, 0.95), weights = c(1, 2, 4),
  prob = c(0.2, 0.5, 0.8)
)
}
\keyword{package}
\keyword{distribution}
