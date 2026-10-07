# weightedBernoulli 0.2.0

- Added `randomized`, `two.sided`, and `draw` arguments to the existing
  `qweighted_bernoulli()` function; no additional function is exported.
- Exact randomization returns the adjacent support cutoffs, their event
  probabilities, the probability of choosing the larger cutoff, and an
  optional independent draw. Log mixing and event probabilities are retained.
- Added the original fair two-sided construction, including the `-1`
  sentinel, the center-atom correction, and exact treatment of zero weights.
- Reused exact-engine probability queries and support information to obtain
  the randomization without full Cartesian-product enumeration.
- Preserved the default exact threshold of 25 active summands and affordable
  structural exceptions. Unresolved larger sums use an explicitly marked
  degenerate conservative cutoff; two-sided fallback preserves support
  membership through a strict greedy subset search.
- `draw = FALSE`, ordinary calls, and degenerate mixtures consume no random
  numbers. Randomized calls always return a mixing-distribution list and
  summary; `details = TRUE` adds the candidate table and diagnostics.
- Fixed direct endpoint CDF comparisons for unequal success probabilities.
  Randomized requests now raise an error if adjacent cutoffs become the
  same representable double after an offset is added.
- Added a mathematical randomization addendum, examples, help documentation,
  and 356 randomized regression checks alongside the 708 ordinary checks.
  The original 77-page note is retained as the version 0.1.0 baseline.

# weightedBernoulli 0.1.0

- Initial R package release with one exported function, `qweighted_bernoulli()`.
- Default exact computation for at most 25 active random summands.
- Exact binomial, repeated weight-probability class, fair rank, and lattice routes.
- Automatic minimum over applicable grid, grouped-count, conditional, and
  concentration bounds within documented work limits.
- Vector probability inputs, logarithmic lower and upper tails, and detailed
  candidate diagnostics.
- Mathematical note, reproducibility material, and installed-package regression
  checks included.
