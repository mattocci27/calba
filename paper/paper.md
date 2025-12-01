---
title: "calba: Efficient neighborhood basal area metrics for trees"
tags:
  - R
  - Spatial Ecology
  - Forest Dynamics
  - Neighborhood Competition
  - Rcpp
authors:
  - name: Masatoshi Katabuchi
    affiliation: 1
affiliations:
  - name: Yunnan Key Laboratory of Forest Ecosystem Stability and Global Change Response, Xishuangbanna Tropical Botanical Garden, Chinese Academy of Sciences, Mengla, Yunnan 666303, China
    index: 1
date: 27 November 2025
bibliography: calba.bib
---

## Summary

`calba` (Calculate Basal Area for Trees) is an R package that compiles neighborhood basal area summaries for every stem in a mapped forest plot.
Built with Rcpp, it computes conspecific and total basal area within configurable radii, applies unweighted or decay-weighted kernels, and flags trees within `r` of the plot boundary via a `safe` edge correction.
The outputs are model-ready covariates for competition, growth, or survival models, and the implementation is tuned to reuse distance information across radii and kernels so large datasets remain tractable.

## Statement of need

Quantifying neighborhood structure is fundamental in spatial ecology and forest dynamics because competition and facilitation operate at local scales defined by distance and species identity [@Song2024; @Hulsmann2024; @Wiegand2025].
Researchers often need to evaluate many alternative hypotheses, such as multiple radii, different decay functions, or focal subsets of species, before fitting statistical models.
Generic spatial packages can compute pairwise distances or intensity functions, but they rarely provide conspecific and total basal-area metrics with built-in edge correction and customizable kernels.
Packages such as `siplab` [@Garcia2014], `forestecology` [@Kim2021], and the ForestGEO package `fgeo` include neighborhood or spatial utilities that support a wide range of tasks.
However, their neighborhood functions are typically designed for broader plot-level workflows and therefore do not focus on producing conspecific and total basal-area summaries with flexible radii, distance-decay kernels, and simple edge treatment in a single, model-ready step.
In contrast, `calba` provides a dedicated high-performance implementation that reuses distance calculations, minimizes memory allocation, and processes radii and decay parameters in a single C++ pass.
These optimizations make `calba` substantially faster in practice, often by orders of magnitude compared with R-level implementations and noticeably faster even than straightforward C++ loops.
`calba` also returns clean, model-ready neighborhood covariates—including conspecific and total basal area, heterospecific components, decay-weighted metrics, and neighbor counts—so that researchers can proceed directly to competition, growth, or survival models without additional data wrangling.
Together, these features allow ecologists to compute comprehensive neighborhood summaries efficiently, even for large forest plots with tens or hundreds of thousands of stems.

## Usage

Install the development version with `remotes::install_github("mattocci27/calba")`.
`neigh_ba()` is the primary interface for computing conspecific and total basal area, neighbor counts, and decay-based summaries in a single call.
The function returns a tidy data frame that can be used directly in regression and mixed-effects models.

```r
library(calba)

set.seed(123)
trees <- data.frame(
  sp = sample(c("oak", "pine", "birch"), 1000, replace = TRUE),
  gx = runif(1000, 0, 100),
  gy = runif(1000, 0, 100),
  ba = runif(1000, 0, 1)
)

neigh_ba(
  sp = trees$sp,
  gx = trees$gx,
  gy = trees$gy,
  ba = trees$ba,
  r = 5,
  mu_values = c(1, 3)
)$summary
```

For analyses that compare multiple radii, `neigh_multi_r()` reuses a single distance matrix and returns neighborhood metrics across user-defined radii.

```r
neigh_multi_r(
  sp = trees$sp,
  gx = trees$gx,
  gy = trees$gy,
  ba = trees$ba,
  r = c(2, 5, 10)
)
```

Advanced users can access the underlying kernels through `ba_simple()`, `ba_decay()`, `count_con()`, and `count_total()`, which provide low-level components for custom workflows.


## Implementation and quality control

All core routines live in `src/calba.cpp`, where distance calculations, species comparisons, and kernel applications are performed in compiled loops.
The R layer exposes user-friendly wrappers (`ba_simple()`, `ba_decay()`, `count_con()`, `count_total()`, `neigh_multi_r()`) with argument validation, consistent unit semantics, and shared edge-foundation logic.
Tests in `tests/testthat/test-calba.R` compare every exported function to brute-force R versions, ensuring the C++ implementation matches expected results across kernels, edge modes, radii, and count helpers.
The package adheres to typical R package checks (`devtools::check()`) and documents the shared library through roxygen tags for reproducible builds.

`calba` includes automated tests using `testthat`, comparing the C++ implementation against brute-force R versions of the same calculations across radii, kernels, and grouping conditions.
The package passes R CMD check with no errors or warnings, and the examples demonstrate correct behavior on a small reproducible dataset.


## Performance and workflow integration

By reusing the same pairwise distance structure for a set of radii or decay parameters, `calba` avoids repeated distance matrix recomputation that would otherwise dominate runtime when exploring model specifications.
This efficiency makes it practical to summarize neighborhood basal area for hundreds or thousands of mapped trees with dozens of kernel/radius combinations.
The outputs can be directly merged into data frames for modeling with GLMMs, GAMs, or Bayesian workflows, and optional helpers simplify counting neighbors or re-running summaries over multiple radii (`neigh_multi_r()`).


## Availability

Install the latest release from GitHub with `remotes::install_github("mattocci27/calba")`.
Source code is hosted at `https://github.com/mattocci27/calba` under the GPL-3 license, and the current version is 0.1.0.
`calba` depends on Rcpp for the compiled kernels and uses testthat for regression tests (`tests/testthat/*`).
Run the bundled tests with `devtools::test()` and validate the package with `devtools::check()` before submission.

## Acknowledgements

## References
