# Validation record: weightedBernoulli 0.2.0

Completed on 7 October 2026 with native **R 4.3.3 on Ubuntu 24.04.3 LTS**.

## Repository-specific refresh

This source archive includes the corrected repository name, package subdirectory,
and relative documentation links for `Arun-Kuchibhotla/Weighted-Bernoulli`.
It was rebuilt and checked after those README changes. The quantile source,
installed help, examples, and both regression suites are identical to the
previously validated version 0.2.0 implementation. The refreshed package check
ran both suites again, and the installed randomized-cutoff example passed
([refresh example log](github-update-smoke.log)).

The existing repository-root workflow at base commit
`910b9686ff35dc6f99aab893f39be00da3aaeef4` was preserved.

## Final package gate

The supplied source archive passed `R CMD build`, staged installation, and
`R CMD check --no-manual` with **0 errors, 0 warnings, and 0 notes**.
Both the original regression suite and the new randomized-quantile suite ran
against the installed package inside that final check.

| Check | Result | Evidence |
|:--|:--|:--|
| Source build and installation | Passed | [build.log](build.log), [install.log](install.log) |
| Final R package check | `Status: OK` | [R-CMD-check.log](R-CMD-check.log) |
| Original ordinary-quantile tests | 708 checks passed | [unified_validation.csv](validation/unified_validation.csv) |
| Independent randomized-quantile tests | 356 checks passed across 18 categories | [randomized_validation.csv](validation/randomized_validation.csv) |
| Documented validation runner | All 1,064 checks passed | [validation.log](validation.log) |
| Installed usage examples | Passed | [usage.log](usage.log), [usage_validation.txt](usage/results/usage_validation.txt) |
| Package and standalone function with `R_DEFAULT_PACKAGES=NULL` | Passed | [minimal_namespace.log](minimal_namespace.log) |
| Exactly one public export, `qweighted_bernoulli` | Passed | [minimal_namespace.log](minimal_namespace.log) |

The `--no-manual` option skips rebuilding an R help-manual PDF. The Rd files,
usage sections, documentation coverage, examples, and both test scripts were
checked. The new six-page mathematical supplement was separately compiled
and visually reviewed on every page.

## What the new tests establish

The randomized tests use independent finite-law enumeration and analytic
references to check actual event probabilities. They cover:

- The original fair two-sided cutoff, its sentinel `-1`, center atoms,
  missing centers, fractional weights, and support gaps.
- One-sided strict-tail randomization with unequal Bernoulli probabilities,
  deterministic terms, endpoint levels, and true support predecessors.
- Equality at CDF jumps, including the unequal-probability endpoint cases
  that exposed avoidable logarithm/exponential roundtrip errors.
- Exact binomial, lattice, repeated-class, and meet-in-the-middle routes,
  including the default length-25 guarantee.
- Logarithmic mixing weights and rejection probabilities beyond ordinary
  floating-point underflow.
- Reproducible sampling, preservation of RNG state with `draw=FALSE`, and
  absence of RNG use in ordinary calls and degenerate mixtures.
- Conservative fallback probabilities evaluated against the actual small
  reference laws, including a winning grouped bound at the required
  two-sided level `alpha/2`.
- A precision error when adding a large deterministic offset collapses two
  adjacent cutoff values to the same representable double.

Probability identities are mathematical statements under independent
randomization. The implementation retains the documented double-precision
and finite-resolution RNG qualifications. The tests establish the stated
numerical comparisons; they are not a general interval-arithmetic certificate.

## Reference results and version history

The original 77-page theoretical note and its numerical reference results
remain the version 0.1.0 baseline. The current randomized API, its proofs,
and its conservative fallback appear in `inst/doc/randomized_quantiles.pdf`.
The current implementation is `R/qweighted_bernoulli.R`.

The version 0.1.0 validation record, including its 42 comparison cases and
261 individual candidate-tail checks, is retained separately under
[docs/validation-v0.1.0](../validation-v0.1.0/README.md). Those historical
comparison counts are not presented as a second execution for this release.
The ordinary implementation's regression suite was rerun for version 0.2.0.

The GitHub Actions workflow is configured for release R on Linux, macOS,
and Windows. The results recorded here are the local native Linux checks;
the hosted jobs execute after the repository files are committed.

## Artifact integrity

Source archive: `weightedBernoulli_0.2.0.tar.gz`.

```text
SHA-256: 3de9e872088cfde2b1a08f799c4f29a4a076b84b98350c6bc485c82045abf9b1
```

Function source, identical in the package archive, package folder, and
standalone download:

```text
SHA-256: b9fd1331cd815af5a00cafd76898bbc0d4e2adba88c8970c8ead7dc40180bac0
```

The validation directory is excluded from R source builds by `.Rbuildignore`.
Adding these records does not change the checked source archive.
