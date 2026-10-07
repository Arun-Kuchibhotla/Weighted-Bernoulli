# Validation record: weightedBernoulli 0.1.0

Validation completed on 7 October 2026 using **native R 4.3.3 on Ubuntu 24.04.3 LTS**, on `x86_64-pc-linux-gnu`.

## Package checks

The final source archive passed:

```sh
R CMD build weightedBernoulli
R CMD INSTALL weightedBernoulli_0.1.0.tar.gz
R CMD check --no-manual weightedBernoulli_0.1.0.tar.gz
```

**Final result: `Status: OK` — 0 errors, 0 warnings, 0 notes.** The checked archive includes the final README and the research-source layout described there.

| Check | Result | Evidence |
|:--|:--|:--|
| Source package build | Passed | [build.log](build.log) |
| Staged installation, byte compilation, and load checks | Passed | [install.log](install.log) |
| `R CMD check --no-manual` | `Status: OK` | [check.log](check.log) |
| Installed namespace with `R_DEFAULT_PACKAGES=NULL` | Passed | [smoke.log](smoke.log) |
| Exactly one public export, `qweighted_bernoulli` | Passed | [regression.csv](regression.csv) |
| `citation()` and installed PDF, TeX, and script discovery | Passed | [smoke.log](smoke.log) |
| Full installed-package regression suite | 708 checks passed across 25 categories | [regression.log](regression.log), [regression.csv](regression.csv) |
| Executable package usage examples | All checks passed | [usage.log](usage.log), [usage_validation.txt](usage_validation.txt) |

The regression suite ran both inside the final package check and through the repository's `tools/run-validation.R` runner. The standalone runner also worked from a directory outside the package source and wrote its report to the requested output directory. Its observed run time was 2.49 seconds; this is one validation run, rather than a portable performance guarantee.

The check used `--no-manual`, so it did not rebuild an R help manual PDF. It did validate the Rd files, usage sections, cross-references, documentation coverage, and help examples. The separately supplied mathematical-note PDF remains included.

The GitHub workflow is configured to check release R on Linux, macOS, and Windows. Those hosted jobs will first execute after the files are committed to the user's GitHub repository. The native results recorded here are for the Linux environment above.

## Installed-package numerical reproduction

The command below was run successfully against the installed package:

```sh
Rscript tools/reproduce.R reproduced-results
```

It completed all 14 settings at three tail levels per setting:

- **42 setting-level comparisons**, comprising 18 exact and 24 conservative returned values.
- **8 settings with independent reference distributions.** Every exact returned value with an available independent reference agreed with that reference quantile.
- **261 individual candidate-tail checks** against reference laws. Every checked strict tail satisfied `bound * (1 + 1e-10) + 1e-15`.
- Every conservative returned value equaled the minimum of the available candidates for that request.
- All comparison tables and both research figures were regenerated successfully.

The numerical tolerance above describes the validation comparison; it does not replace the mathematical real-arithmetic statements or certify every floating-point operation.

The settings include lengths 25 and 26 around the default exact threshold, equal weights, fair rank weights, a dominant weight, a heterogeneous lattice, two repeated real-weight classes, square-root and harmonic weights, nearly equal weights at lengths 100 and 1,000, unequal success probabilities, and rare success probabilities.

Evidence from this native run is included as [comparisons.csv](comparisons.csv), [all_bounds.csv](all_bounds.csv), [comparison_validation.txt](comparison_validation.txt), and [native_sessionInfo.txt](native_sessionInfo.txt).

The original mathematical note and its adjacent reference results under `inst/research/` retain their documented R/WebAssembly environment and timings. The files in this directory record the separate native package run. No timing claim in the note has been silently replaced with a timing from this machine.

## Source integrity

The packaged `R/qweighted_bernoulli.R` is byte-for-byte identical to the validated standalone implementation supplied before packaging. Its SHA-256 is:

```text
e42e441541f3f55f539ffc62b009022a975729933d7bb330392d9484031295d2
```

The final, checked `weightedBernoulli_0.1.0.tar.gz` has SHA-256:

```text
01daf3183d3c21fc67afa07615773d9495c6977470b028f1ddb760bf1aa1dc59
```

All computational helpers remain local to the single exported function. Packaging adds metadata, a namespace import for `utils::tail`, installed help, citation metadata, package-aware example and reproduction entry points, tests, and GitHub automation.

The `docs/` directory is included in the GitHub folder and excluded from the R source build by `.Rbuildignore`. Adding this validation record therefore does not change the checked source archive.
