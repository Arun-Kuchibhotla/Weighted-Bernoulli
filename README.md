# weightedBernoulli

**Version 0.2.0 adds randomized quantiles, including the original fair two-sided cutoffs `c_plus`, `c_minus`, and `tau`, to the same single R function.** Exact quantiles and the automatic minimum of applicable conservative bounds remain available with the original defaults.

For independent random variables $\xi_i \sim \operatorname{Bernoulli}(\pi_i)$, the package computes quantiles of

$$
G = \sum_{i=1}^{B} w_i \xi_i,
\qquad w_i \geq 0.
$$

The only exported function is `qweighted_bernoulli()`. It accepts a common success probability or a separate probability for each weight, and handles several quantile levels in one call. Use `randomized = TRUE` to obtain the mixing distribution and an optional independent draw; use `two.sided = TRUE` as well for the original symmetric inclusive-tail cutoff.

- **At most 25 active summands:** compute the exact finite-distribution quantile by default.
- **Larger inputs with affordable exact structure:** still compute the exact quantile. Examples include equal weights, repeated weight–probability classes, fair rank weights, and verified lattices.
- **Other larger inputs:** compute every applicable conservative candidate within the documented work limits and return the **minimum** of those candidates.

An active summand has a positive weight and a success probability strictly between zero and one. Zero terms are removed; deterministic terms contribute an offset. “Exact” refers to the finite-distribution algorithm evaluated in double precision.

## Installation

Requires **R 4.1.0 or later**. Runtime dependencies are limited to the standard `stats` and `utils` packages supplied with R. The package contains no compiled code.

### From GitHub

Install from [Arun-Kuchibhotla/Weighted-Bernoulli](https://github.com/Arun-Kuchibhotla/Weighted-Bernoulli), specifying the package subdirectory:

```r
install.packages("remotes")  # Once, if needed.

remotes::install_github(
  "Arun-Kuchibhotla/Weighted-Bernoulli",
  subdir = "weightedBernoulli"
)

library(weightedBernoulli)
```

If an older version is already loaded, restart R before reinstalling. The [`remotes` documentation](https://remotes.r-lib.org/reference/install_github.html) describes branch and release selection.

### From the downloaded folder

On the repository's GitHub page, choose **Code > Download ZIP**, then extract the download. For the `main` branch, the outer folder is named `Weighted-Bernoulli-main`. Open a terminal in that folder—the repository root, containing the `weightedBernoulli` subfolder—and run:

```sh
R CMD INSTALL weightedBernoulli
```

Then, in R:

```r
library(weightedBernoulli)
```

Alternatively, set R's working directory to the repository root and install the built source archive directly:

```r
install.packages(
  "weightedBernoulli_0.2.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

Supply the full path if the archive is outside your current working directory. Installing this package does not require a C/C++ compiler or LaTeX.

## Randomized cutoffs

### The original fair two-sided construction

For fair Bernoulli variables, let $N = \sum_i w_i$ and define

$$
P_h(c) = \Pr\left(\left|G-N/2\right| \geq N/2-c\right),
\qquad P_h(-1)=0.
$$

The exact branch returns the first support point `c_plus` at which $P_h(c)\geq\alpha$, its support predecessor `c_minus` (using the sentinel `-1` below zero), and the mixing probability

$$
\tau = \frac{\alpha-P_h(c_-)}{P_h(c_+)-P_h(c_-)}.
$$

Here `tau` is the probability of choosing the **larger** cutoff `c_plus`. First compute the distribution without drawing from it:

```r
library(weightedBernoulli)

cut <- qweighted_bernoulli(
  p = 0.05,
  weights = 1:20,
  lower.tail = FALSE,
  randomized = TRUE,
  two.sided = TRUE,
  draw = FALSE
)

cut$c_minus
# [1] 52
cut$c_plus
# [1] 53
round(cut$tau, 10)
# [1] 0.3297297297
cut$kappa
# NULL

cut$tau * cut$P_plus + (1 - cut$tau) * cut$P_minus
# [1] 0.05
cut$exact
# [1] TRUE
```

To generate the cutoff as well, use the default `draw = TRUE`:

```r
set.seed(20261007)
sampled <- qweighted_bernoulli(
  p = 0.05,
  weights = 1:20,
  lower.tail = FALSE,
  randomized = TRUE,
  two.sided = TRUE
)

sampled$kappa  # One draw: either 52 or 53.
sampled$summary
```

For an independent fresh uniform randomizer, the exact mixture satisfies

$$
\Pr\left(\left|G-N/2\right|\geq N/2-\kappa\right)=\alpha.
$$

The probability includes both the Bernoulli sum and the independent cutoff randomization. Do not reuse a tie-breaking uniform that participates in the statistic. Two-sided mode requires `randomized = TRUE`, `prob = 0.5` on every positive weight, and a tail level strictly between zero and one. Supplying `p = 0.95` with `lower.tail = TRUE` specifies the same level as `p = 0.05` with `lower.tail = FALSE`.

**An atom at the center is counted once.** For $c<N/2$, symmetry gives $P_h(c)=2\Pr(G\leq c)$, whereas $P_h(N/2)=1$. For example:

```r
center <- qweighted_bernoulli(
  p = 0.75, weights = c(1, 1), lower.tail = FALSE,
  randomized = TRUE, two.sided = TRUE, draw = FALSE
)
center[c("c_minus", "c_plus", "tau", "P_minus", "P_plus")]
# c_minus = 0, c_plus = 1, tau = 0.5
# P_minus = 0.5, P_plus = 1
```

### One-sided randomization, including unequal probabilities

Without `two.sided = TRUE`, the mixture calibrates the **strict upper tail**. In the exact branch, `c_plus` is the ordinary $(1-\alpha)$ quantile and `c_minus` is its support predecessor, with `-Inf` below the minimum. Here the probability fields are $P_\pm=\Pr(G>c_\pm)$ and

$$
\tau=\frac{P_- - \alpha}{P_- - P_+},
\qquad \Pr(G>\kappa)=\alpha.
$$

```r
one <- qweighted_bernoulli(
  p = 0.05, weights = c(1, 2, 4),
  prob = c(0.2, 0.5, 0.8), lower.tail = FALSE,
  randomized = TRUE, draw = FALSE
)
one[c("c_minus", "c_plus", "tau", "P_minus", "P_plus")]
# c_minus = 6, c_plus = 7, tau = 0.375
# P_minus = 0.08, P_plus = 0
```

Although the two formulas have different orientations, `tau` always selects `c_plus`. The direction changes because increasing a one-sided upper threshold decreases its strict-tail probability, whereas increasing the original two-sided cutoff increases its event probability. At ordinary lower-tail level zero, one-sided randomization returns `-Inf`; at level one, it returns the maximum support point. The ordinary nonrandomized function retains its support-endpoint convention.

### When exact randomization is too expensive

The function preserves the exact-short-vector and structural-exception policy. If it must use conservative bounds, it reports a **degenerate conservative mixture**: `c_plus == c_minus`, `tau = 1`, and `exact = FALSE`. In those rows, `P_plus` and `P_minus` are **upper bounds**, and the weighted probability calculation certifies `<= alpha`. They do not identify the exact adjacent bracket or claim exact nominal attainment.

For one-sided mode, this uses the minimum available deterministic quantile bound. For two-sided mode, the function first obtains the best available strict upper-tail bound `q` at level `alpha/2`, then finds a feasible subset sum strictly below `N - q`; the sentinel `-1` is used when necessary. Ascending and descending greedy subset searches preserve actual support membership without enumerating the full distribution. The strict inequality prevents an uncontrolled boundary atom from entering the inclusive two-sided event. Candidate selection takes place before sampling.

Randomized output always includes its mixing distribution and summary. With `draw = FALSE`, both `kappa` and its alias `quantile` are `NULL`, and the RNG state is unchanged. Ordinary calls and degenerate mixtures also use no RNG. `details = TRUE` adds the candidate table and diagnostics. The [randomization addendum](weightedBernoulli/inst/doc/randomized_quantiles.pdf) provides the definitions, proofs, computation details, and conservative guarantees.

## Deterministic quantiles: quick start

### Exact quantiles with unequal success probabilities

```r
library(weightedBernoulli)

qweighted_bernoulli(
  p = c(0.50, 0.90, 0.95),
  weights = c(1, 2, 4),
  prob = c(0.2, 0.5, 0.8)
)
# [1] 4 6 7
```

Here `p` contains the requested quantile levels and `prob` contains the Bernoulli success probabilities.

### Automatic exact calculation for short vectors

```r
fit <- qweighted_bernoulli(
  p = 0.95,
  weights = sqrt(1:25),
  details = TRUE
)

fit$quantile
# Approximately 57.66373

fit$summary[, c("quantile", "exact", "method")]
# exact is TRUE; the general exact calculation uses meet-in-the-middle.
```

The default is `max_length = 25`; no method argument is needed. For a general length-25 problem, meet-in-the-middle uses half-distributions with 4,096 and 8,192 states. If a required exact calculation cannot resolve a CDF jump within its numerical or work limits, it raises an error.

### Exact structure above the length threshold

```r
# One large weight and 39 equal weights reduce to two binomial counts.
fit <- qweighted_bernoulli(
  p = 0.05,
  weights = c(60, rep(1, 39)),
  lower.tail = FALSE,
  details = TRUE
)

fit$quantile
# [1] 83
fit$summary$exact
# [1] TRUE

# A thousand equal weights reduce to one binomial distribution.
qweighted_bernoulli(0.95, weights = rep(1, 1000))
# [1] 526

# Fair rank weights use their exact signed-rank distribution.
qweighted_bernoulli(0.95, weights = 1:100)
# [1] 3004
```

The threshold concerns generic exact subset-sum calculations. Affordable structural calculations are considered at any length, including repeated classes with different success probabilities across classes.

### Compare the bounds for a longer vector

```r
fit <- qweighted_bernoulli(
  p = c(0.05, 0.01, 1e-6),
  weights = sqrt(1:100),
  lower.tail = FALSE,
  details = TRUE
)

round(fit$quantile, 4)
# [1] 394.3164 418.2969 501.0234

fit$summary

# All available candidates for the first requested level, best first.
tab <- subset(fit$bounds, request == 1 & available)
tab <- tab[order(tab$quantile), ]
tab[, c("method", "quantile", "selected")]

stopifnot(fit$quantile[1] == min(tab$quantile))

# Reasons that individual candidates were unavailable.
subset(
  fit$bounds,
  !available,
  select = c("request", "method", "reason")
)
```

This example returns conservative upper bounds. Multiple levels in one call share preparation of distributions and bounds. `details = TRUE` changes the returned diagnostics; it does not change which candidates are evaluated.

## Quantile and tail conventions

For $0 < u < 1$, the exact quantile is

$$
q_G(u) = \inf\{x : \Pr(G \leq x) \geq u\}.
$$

A conservative result $\widehat q(u)$ is an upper bound for $q_G(u)$. With `lower.tail = FALSE` and an upper-tail level $\alpha$, the target is the $(1-\alpha)$ quantile, with the guarantee

$$
\Pr\{G > \widehat q(1-\alpha)\} \leq \alpha.
$$

The tail event uses **`>`**. An inclusive event such as `G >= q` also includes any probability mass at `q`, so a test with an inclusive rejection rule must account for that atom. A conservative threshold need not be a support point. At quantile levels zero and one, the ordinary nonrandomized function returns the lower and upper support endpoints, respectively.

For very small upper-tail probabilities, supply the tail directly. Using `log.p = TRUE` also avoids ordinary-scale underflow:

```r
qweighted_bernoulli(
  p = log(1e-30),
  weights = sqrt(1:26),
  lower.tail = FALSE,
  log.p = TRUE
)
```

Forming `1 - 1e-30` first would lose the tail level in double precision.

### Signed weights

The function accepts nonnegative weights. A signed-weight problem reduces exactly to this interface by replacing each negative-weight Bernoulli variable with its complement and adding an offset:

```r
weights <- c(-3, 1, 2)
prob <- c(0.2, 0.5, 0.8)
negative <- weights < 0

offset <- sum(weights[negative])
prob[negative] <- 1 - prob[negative]

offset + qweighted_bernoulli(
  p = 0.95,
  weights = abs(weights),
  prob = prob
)
# [1] 3
```

## Interface

```r
qweighted_bernoulli(
  p,
  weights,
  prob = 0.5,
  max_length = 25,
  lower.tail = TRUE,
  log.p = FALSE,
  details = FALSE,
  control = list(),
  randomized = FALSE,
  two.sided = FALSE,
  draw = TRUE
)
```

| Argument | Meaning |
|:--|:--|
| `p` | One or more quantile probabilities; upper-tail levels if `lower.tail = FALSE`; logarithms if `log.p = TRUE`. Names are preserved. |
| `weights` | Finite, nonnegative weights with a finite sum. Empty vectors and zero weights are allowed. |
| `prob` | A single Bernoulli success probability, or one probability per weight; values in `[0, 1]`. |
| `max_length` | Maximum active length for generic exact calculation; default `25`, allowed integers `0` through `52`. Structural exact calculations are checked at any length. |
| `lower.tail` | Use lower-tail quantile probabilities when `TRUE`; use upper-tail levels when `FALSE`. |
| `log.p` | Interpret `p` on the natural logarithmic scale. |
| `details` | Add the candidate table and diagnostics; ordinary calls also gain a summary. Randomized calls always include their mixing distribution and summary. |
| `control` | Named overrides of computational work limits; see below and the installed help page. |
| `randomized` | If `TRUE`, return a mixing distribution and, by default, a draw. The default `FALSE` retains the ordinary interface. |
| `two.sided` | If `TRUE`, compute the original inclusive two-sided cutoff for fair Bernoulli sums. Requires `randomized = TRUE` and `0 < alpha < 1`. |
| `draw` | In randomized mode, generate `kappa` when `TRUE`; return the distribution without using the RNG when `FALSE`. |

For ordinary calls (`randomized = FALSE`), `details = FALSE` returns a numeric vector. With `details = TRUE`, the result has these components:

| Component | Contents |
|:--|:--|
| `quantile` | The same numeric quantile vector. |
| `summary` | One row per level: the selected quantile, whether it is exact, the method, strict-tail bound fields, and a grid allowance where applicable. |
| `bounds` | Candidate values, availability, selection, reasons for unavailable candidates, partitions, and state counts. `NULL` if every request was resolved exactly. |
| `diagnostics` | Input and active lengths, completed work controls, exact attempts, and partition information. |

The summary's `tail_bound` and `log_tail_bound` are upper bounds at the returned threshold; they do not report a newly measured achieved tail probability. An exact interior result usually records the requested upper-tail level. `grid_error` is zero for exact answers, an absolute allowance for selected grid or grouped-envelope candidates, and `NA` for other selected bounds. The help page documents every column.

For randomized calls, the result is always a list:

| Component | Contents |
|:--|:--|
| `c_plus`, `c_minus` | Exact adjacent cutoffs, or equal endpoints for a degenerate conservative fallback. |
| `tau` | Probability of choosing `c_plus`. |
| `log_tau`, `log_one_minus_tau` | Log mixing probabilities, retained when ordinary-scale values round to zero or one. |
| `P_plus`, `P_minus` | Exact event probabilities when `exact = TRUE`; upper bounds when `exact = FALSE`. The event depends on `two.sided`, as defined above. |
| `log_P_plus`, `log_P_minus` | Corresponding natural logarithms. |
| `kappa`, `quantile` | The same generated cutoff vector; both are `NULL` if `draw = FALSE`. |
| `exact` | Whether the exact-law mixing distribution was resolved, one flag per request. |
| `alpha`, `log_alpha`, `two.sided` | Requested tail levels and the selected event convention. |
| `summary` | One row per request, including cutoffs, `tau`, probability fields, `exact`, `degenerate`, `probability_type`, and `method`. |
| `bounds`, `diagnostics` | Added only with `details = TRUE`; describe the underlying exact/bound calculations. |

For conservative two-sided rows, `bounds` contains the underlying upper quantile candidates at level `alpha/2`, before reflection and the support search. Their `quantile` column therefore differs from the final two-sided `c_plus`.

## Methods and computational limits

The conservative collection includes support, symmetry, upward-grid, Chernoff, Hoeffding, Cantelli, Bernstein, and binomial $P_2$ bounds. Fair Bernoulli sums additionally admit Bentkus–Dzindzalieta and Pinelis Gaussian comparisons and Montgomery-Smith $K$-functional refinements. Count aggregation adds grouped envelopes and conditional support, moment, Serfling, and finite-grid MGF calculations. Proofs and comparisons appear in the [mathematical note](weightedBernoulli/inst/doc/weighted_bernoulli_quantiles.pdf).

Count-based conservative refinements currently require a common success probability among the active summands. Other applicable candidates remain available with unequal probabilities. The function records omitted or unsuccessful candidates and returns the minimum over the implemented candidates completed within the work limits. The search covers a documented collection of partitions and MGF parameters; it does not optimize over every possible partition or probability inequality.

The 11 recognized controls and their defaults are:

| Control | Default | Main use |
|:--|--:|:--|
| `max_states` | `2^18` | Lattice/grid arrays and one exact half-distribution. |
| `max_dp_work` | `2e7` | Lattice/grid dynamic-programming updates. |
| `max_mitm_work` | `2e7` | Exact half-distribution query work. |
| `max_iterations` | `80` | Scalar searches and inversions. |
| `chernoff_iterations` | `60` | Full-MGF Chernoff optimization per distinct level. |
| `max_p2_B` | `10000` | Active count for the optional binomial $P_2$ candidate. |
| `max_p2_work` | `200000` | Binomial $P_2$ preparation and scan work. |
| `max_group_states` | `65536` | Non-pivot states in grouped-count calculations. |
| `max_count_states` | `100000` | Joint count states for conditional calculations. |
| `max_mgf_work` | `2e7` | Conditional MGF preparation and table work per partition. |
| `mgf_points` | `37` | Total conditional MGF parameters, including zero. |

For example:

```r
fit <- qweighted_bernoulli(
  0.95,
  weights = sqrt(1:100),
  details = TRUE,
  control = list(max_states = 2^19, max_dp_work = 4e7)
)
```

These controls are operation and state-count proxies. Changing them can change the automatically constructed grids and partitions, so a larger budget does not promise a smaller returned bound. For a required exact calculation with at most 25 active terms, the implementation raises the exact storage and search floors and bypasses the exact query-work cap. Lowering conservative budgets therefore preserves the default short-vector exact policy. Larger user-selected values of `max_length` remain subject to the specified exact work limits.

### Numerical scope

The mathematical guarantees are proved in real arithmetic. The implementation uses R double precision, direct small-tail calculations, logarithmic probabilities where needed, directional rounding safeguards, and feasible upper search brackets. This is an exact-distribution or conservative-bound algorithm according to the selected route, with the usual floating-point qualification; it is not symbolic or certified interval arithmetic. Extremely close subset sums and requests at numerically indistinguishable CDF jumps require this qualification. A randomized request raises an error if the adjacent cutoffs collapse to the same representable number. Version 0.2.0 also preserves direct endpoint comparisons for unequal Bernoulli probabilities at CDF jumps.

Draws use the currently selected R generator through `runif()`. Its finite resolution introduces the usual small discretization error relative to an ideal continuous uniform; exceptionally small mixing probabilities may fall below that resolution even when their logarithms are accurately recorded. Use `draw = FALSE` to retain the distribution for a different sampling implementation when that distinction matters. See the official R documentation for [random number generation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html) and [`runif()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Uniform.html).

## Documentation and reproducibility

- [Randomized-quantile addendum (PDF)](weightedBernoulli/inst/doc/randomized_quantiles.pdf): the version 0.2.0 randomized interface, exact identities, center-atom correction, and conservative fallback, with worked examples.
- [Addendum source](weightedBernoulli/inst/research/randomized_quantiles.tex): self-contained LaTeX source.
- [Main mathematical note (PDF)](weightedBernoulli/inst/doc/weighted_bernoulli_quantiles.pdf): the 77-page version 0.1.0 baseline, with derivations, proofs, citations, interlaced R code, numerical comparisons, and the original randomized-cutoff discussion.
- [Main-note source](weightedBernoulli/inst/research/weighted_bernoulli_quantiles.tex): the baseline note's source, with supporting figures and tables in the adjacent folders.
- [Executable usage examples](weightedBernoulli/inst/examples/usage_examples.R): package-based examples with checks of their results.
- [Comparison script](weightedBernoulli/inst/reproduce/reproduce.R): regenerates the comparison tables and figures using the installed package.
- [Ordinary regression suite](weightedBernoulli/tests/regression.R): the existing 708 checks covering the API, exact laws, conservative candidates, dispatch, small tails, and numerical edge cases.
- [Randomized regression suite](weightedBernoulli/tests/randomized.R): 356 additional checks for mixing probabilities, inclusive and strict event conventions, exact-engine metadata, endpoints, RNG behavior, and conservative fallback.
- [Validation record](weightedBernoulli/docs/validation/README.md): the package build, installation, and test results for this distribution.

The main note was written around the version 0.1.0 standalone function. Its printed code appendix is retained as a baseline; **it is not the current implementation**. Use [the package source](weightedBernoulli/R/qweighted_bernoulli.R) for the current function, and the addendum for the randomized API and guarantees. The original ordinary calls remain supported. When following the baseline note with the package installed, replace its initial `source("qweighted_bernoulli.R")` call with `library(weightedBernoulli)`. The packaged comparison and example scripts already make this change.

The paths and terminal commands in the following table are relative to the repository root:

| Filename used in the standalone note | Repository location or command |
|:--|:--|
| `qweighted_bernoulli.R` | `weightedBernoulli/R/qweighted_bernoulli.R`; load with `library(weightedBernoulli)`. |
| `usage_examples.R` | `weightedBernoulli/inst/examples/usage_examples.R`. |
| `reproduce.R` | `Rscript weightedBernoulli/tools/reproduce.R [output_directory]`. |
| `verify_quantiles.R` | `weightedBernoulli/tests/regression.R`, or `Rscript weightedBernoulli/tools/run-validation.R [output_directory]`. |

From R:

```r
help("qweighted_bernoulli", package = "weightedBernoulli")
citation("weightedBernoulli")

system.file(
  "doc", "randomized_quantiles.pdf",
  package = "weightedBernoulli",
  mustWork = TRUE
)
```

To run the installed examples in a writable output directory:

```r
example_script <- system.file(
  "examples", "usage_examples.R",
  package = "weightedBernoulli",
  mustWork = TRUE
)
example_output <- file.path(tempdir(), "weightedBernoulli-examples")
dir.create(example_output, recursive = TRUE, showWarnings = FALSE)
old_directory <- setwd(example_output)
tryCatch(source(example_script), finally = setwd(old_directory))
```

After installing the package, run these commands from the repository root to reproduce the research comparisons or write the regression report:

```sh
Rscript weightedBernoulli/tools/reproduce.R reproduced-results
Rscript weightedBernoulli/tools/run-validation.R validation-results
```

The first command creates tables, plots, and its own `sessionInfo.txt`; the second creates `unified_validation.csv` and `randomized_validation.csv` for the 708 ordinary and 356 randomized checks (1,064 in total). The shipped note and its embedded timings retain the original documented runtime. Reproduction on another R version or machine writes that environment's results into a separate directory.

To rebuild either document, change from the repository root to `weightedBernoulli/inst/research` and run `pdflatex` twice on `weighted_bernoulli_quantiles.tex` or `randomized_quantiles.tex`; LaTeX is needed only for that optional documentation step. `R CMD check` runs both the ordinary and randomized test scripts.

## Repository layout and development

The repository is `Arun-Kuchibhotla/Weighted-Bernoulli`; the R package lives in its `weightedBernoulli/` subdirectory. All terminal commands below assume the repository root as the working directory.

1. Clone the repository or download and extract it using **Code > Download ZIP**.
2. Keep `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, and the other package files inside `weightedBernoulli/`. When updating from a supplied package folder, merge it into the existing folder of that name; the resulting source path is `weightedBernoulli/R/qweighted_bernoulli.R`.
3. Keep the GitHub Actions workflow at `.github/workflows/R-CMD-check.yaml` in the repository root. Its dependency and package-check actions set `working-directory: weightedBernoulli` to find the package. Include the package's `.Rbuildignore` and other repository configuration files in commits.
4. Commit changes to `main` or `master`, or open a pull request. The included workflow runs `R CMD check` with release R on Linux, macOS, and Windows.
5. Install the updated package with the GitHub command above, keeping `subdir = "weightedBernoulli"`.

For a local package check, run from the repository root:

```sh
R CMD build weightedBernoulli
R CMD check --no-manual weightedBernoulli_0.2.0.tar.gz
```

The package uses handwritten R help files and base-R tests, so contributing does not require `roxygen2` or `testthat`. Package conventions follow the official [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html) manual; the workflow follows the [r-lib actions examples](https://github.com/r-lib/actions/tree/v2/examples).

All paths in this table are relative to the repository root:

| Path | Purpose |
|:--|:--|
| `weightedBernoulli/R/qweighted_bernoulli.R` | The one public function, with all computational helpers local to it. |
| `weightedBernoulli/DESCRIPTION`, `weightedBernoulli/NAMESPACE` | Package metadata, runtime imports, and the single export. |
| `weightedBernoulli/man/` | Installed function and package help. |
| `weightedBernoulli/inst/doc/` | Installed baseline mathematical note and randomized-quantile addendum (PDFs). |
| `weightedBernoulli/inst/research/` | LaTeX source, figures, tables, and reference results. |
| `weightedBernoulli/inst/examples/`, `weightedBernoulli/inst/reproduce/` | Installed examples and comparison script. |
| `weightedBernoulli/tests/`, `weightedBernoulli/tools/` | Regression checks and repository command-line runners. |
| `.github/workflows/` | Checks on three operating systems, configured for the package subdirectory. |
| `weightedBernoulli/docs/validation/` | Recorded validation of the supplied package. |

## References

The note contains the full bibliography and identifies the assumptions and constants used by each bound. Selected primary references are:

- Bentkus, V. (2004). *On Hoeffding's inequalities*. The Annals of Probability 32(2), 1650–1673. [doi:10.1214/009117904000000360](https://doi.org/10.1214/009117904000000360).
- Bentkus, V. K., and Dzindzalieta, D. (2015). *A tight Gaussian bound for weighted sums of Rademacher random variables*. Bernoulli 21(2), 1231–1237. [doi:10.3150/14-BEJ603](https://doi.org/10.3150/14-BEJ603).
- Pinelis, I. (2012). *An asymptotically Gaussian bound on the Rademacher tails*. Electronic Journal of Probability 17, article 35. [doi:10.1214/EJP.v17-2026](https://doi.org/10.1214/EJP.v17-2026).
- Montgomery-Smith, S. J. (1990). *The distribution of Rademacher sums*. Proceedings of the American Mathematical Society 109(2), 517–522. [doi:10.1090/S0002-9939-1990-1013975-0](https://doi.org/10.1090/S0002-9939-1990-1013975-0).
- Bardenet, R., and Maillard, O.-A. (2015). *Concentration inequalities for sampling without replacement*. Bernoulli 21(3), 1361–1385. [doi:10.3150/14-BEJ605](https://doi.org/10.3150/14-BEJ605).

## License

Copyright 2026 Arun Kumar Kuchibhotla. Distributed under the [MIT license](weightedBernoulli/LICENSE.md).
