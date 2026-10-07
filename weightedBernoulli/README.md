# weightedBernoulli

**Exact quantiles and conservative upper quantile bounds for sums of independent weighted Bernoulli random variables, through one R function.**

For independent random variables $\xi_i \sim \mathrm{Bernoulli}(\pi_i)$, the package computes quantiles of

$$
G = \sum_{i=1}^{B} w_i \xi_i,
\qquad w_i \geq 0.
$$

The only exported function is `qweighted_bernoulli()`. It accepts a common success probability or a separate probability for each weight, and handles several quantile levels in one call.

- **At most 25 active summands:** compute the exact finite-distribution quantile by default.
- **Larger inputs with affordable exact structure:** still compute the exact quantile. Examples include equal weights, repeated weight–probability classes, fair rank weights, and verified lattices.
- **Other larger inputs:** compute every applicable conservative candidate within the documented work limits and return the **minimum** of those candidates.

An active summand has a positive weight and a success probability strictly between zero and one. Zero terms are removed; deterministic terms contribute an offset. “Exact” refers to the finite-distribution algorithm evaluated in double precision.

## Installation

Requires **R 4.1.0 or later**. Runtime dependencies are limited to the standard `stats` and `utils` packages supplied with R. The package contains no compiled code.

### From GitHub

Install directly from GitHub::

```r
install.packages("remotes")  # Once, if needed.

remotes::install_github(
  "Arun-Kuchibhotla/Weighted-Bernoulli",
  subdir = "weightedBernoulli"
)

library(weightedBernoulli)
```

If you upload it under another GitHub account or repository name, change the `"owner/repository"` string accordingly. The GitHub command becomes usable once the repository has been created and the package files committed. The [`remotes` documentation](https://remotes.r-lib.org/reference/install_github.html) describes branch and release selection.

### From the downloaded folder

Extract `weightedBernoulli-github.zip`. In a terminal, change to the directory containing the extracted `weightedBernoulli` folder and run:

```sh
R CMD INSTALL weightedBernoulli
```

Then, in R:

```r
library(weightedBernoulli)
```

Alternatively, install the accompanying built source archive directly from R:

```r
install.packages(
  "weightedBernoulli_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

Supply the full path if the archive is outside your current working directory. Installing this package does not require a C/C++ compiler or LaTeX.

## Quick start

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

The tail event uses **`>`**. An inclusive event such as `G >= q` also includes any probability mass at `q`, so a test with an inclusive rejection rule must account for that atom. A conservative threshold need not be a support point. At quantile levels zero and one, the function returns the lower and upper support endpoints, respectively.

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
  control = list()
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
| `details` | Return the quantiles together with the summary, candidate table, and diagnostics. |
| `control` | Named overrides of computational work limits; see below and the installed help page. |

With `details = FALSE`, the result is a numeric vector. With `details = TRUE`, the result has these components:

| Component | Contents |
|:--|:--|
| `quantile` | The same numeric quantile vector. |
| `summary` | One row per level: the selected quantile, whether it is exact, the method, strict-tail bound fields, and a grid allowance where applicable. |
| `bounds` | Candidate values, availability, selection, reasons for unavailable candidates, partitions, and state counts. `NULL` if every request was resolved exactly. |
| `diagnostics` | Input and active lengths, completed work controls, exact attempts, and partition information. |

The summary's `tail_bound` and `log_tail_bound` are upper bounds at the returned threshold; they do not report a newly measured achieved tail probability. An exact interior result usually records the requested upper-tail level. `grid_error` is zero for exact answers, an absolute allowance for selected grid or grouped-envelope candidates, and `NA` for other selected bounds. The help page documents every column.

## Methods and computational limits

The conservative collection includes support, symmetry, upward-grid, Chernoff, Hoeffding, Cantelli, Bernstein, and binomial $P_2$ bounds. Fair Bernoulli sums additionally admit Bentkus–Dzindzalieta and Pinelis Gaussian comparisons and Montgomery-Smith $K$-functional refinements. Count aggregation adds grouped envelopes and conditional support, moment, Serfling, and finite-grid MGF calculations. Proofs and comparisons appear in the [mathematical note](inst/doc/weighted_bernoulli_quantiles.pdf).

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

The mathematical guarantees are proved in real arithmetic. The implementation uses R double precision, direct small-tail calculations, logarithmic probabilities where needed, directional rounding safeguards, and feasible upper search brackets. This is an exact-distribution or conservative-bound algorithm according to the selected route, with the usual floating-point qualification; it is not symbolic or certified interval arithmetic. Extremely close subset sums and requests at numerically indistinguishable CDF jumps require this qualification. See the note and help page for details.

## Documentation and reproducibility

- [Mathematical note (PDF)](inst/doc/weighted_bernoulli_quantiles.pdf): a self-contained, 77-page treatment with derivations, proofs, citations, interlaced R code, numerical comparisons, and the original randomized-cutoff discussion.
- [LaTeX source](inst/research/weighted_bernoulli_quantiles.tex): the note's source, with supporting figures and tables in the adjacent folders.
- [Executable usage examples](inst/examples/usage_examples.R): package-based examples with checks of their results.
- [Comparison script](inst/reproduce/reproduce.R): regenerates the comparison tables and figures using the installed package.
- [Regression suite](tests/regression.R): 708 checks covering the API, exact laws, conservative candidates, dispatch, small tails, and numerical edge cases.
- [Validation record](docs/validation/README.md): the package build, installation, and test results for this distribution.

The note was written around the standalone function. When following it with the package installed, replace its initial `source("qweighted_bernoulli.R")` call with `library(weightedBernoulli)`. The function interface and calculations are the same. The packaged comparison and example scripts already make this change.

| Filename used in the standalone note | Package location or command |
|:--|:--|
| `qweighted_bernoulli.R` | `R/qweighted_bernoulli.R`; load with `library(weightedBernoulli)`. |
| `usage_examples.R` | `inst/examples/usage_examples.R`. |
| `reproduce.R` | `Rscript tools/reproduce.R [output_directory]`. |
| `verify_quantiles.R` | `tests/regression.R`, or `Rscript tools/run-validation.R [output_directory]`. |

From R:

```r
help("qweighted_bernoulli", package = "weightedBernoulli")
citation("weightedBernoulli")

system.file(
  "doc", "weighted_bernoulli_quantiles.pdf",
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
Rscript tools/reproduce.R reproduced-results
Rscript tools/run-validation.R validation-results
```

The first command creates tables, plots, and its own `sessionInfo.txt`; the second creates `unified_validation.csv`. The shipped note and its embedded timings retain the original documented runtime. Reproduction on another R version or machine writes that environment's results into a separate directory.

To rebuild the note, change to `inst/research` and run `pdflatex` twice on `weighted_bernoulli_quantiles.tex`; LaTeX is needed only for that optional documentation step.


## References

The note contains the full bibliography and identifies the assumptions and constants used by each bound. Selected primary references are:

- Bentkus, V. (2004). *On Hoeffding's inequalities*. The Annals of Probability 32(2), 1650–1673. [doi:10.1214/009117904000000360](https://doi.org/10.1214/009117904000000360).
- Bentkus, V. K., and Dzindzalieta, D. (2015). *A tight Gaussian bound for weighted sums of Rademacher random variables*. Bernoulli 21(2), 1231–1237. [doi:10.3150/14-BEJ603](https://doi.org/10.3150/14-BEJ603).
- Pinelis, I. (2012). *An asymptotically Gaussian bound on the Rademacher tails*. Electronic Journal of Probability 17, article 35. [doi:10.1214/EJP.v17-2026](https://doi.org/10.1214/EJP.v17-2026).
- Montgomery-Smith, S. J. (1990). *The distribution of Rademacher sums*. Proceedings of the American Mathematical Society 109(2), 517–522. [doi:10.1090/S0002-9939-1990-1013975-0](https://doi.org/10.1090/S0002-9939-1990-1013975-0).
- Bardenet, R., and Maillard, O.-A. (2015). *Concentration inequalities for sampling without replacement*. Bernoulli 21(3), 1361–1385. [doi:10.3150/14-BEJ605](https://doi.org/10.3150/14-BEJ605).

## License

Copyright 2026 Arun Kumar Kuchibhotla. Distributed under the [MIT license](LICENSE.md).
