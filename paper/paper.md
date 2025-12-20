---
title: "calba: An R package for efficient neighborhood basal area metrics for trees"
tags:
  - R
  - Spatial Ecology
  - Forest Dynamics
  - Neighborhood Competition
  - Rcpp
authors:
  - name: Masatoshi Katabuchi
    affiliation: "1,2"
affiliations:
  - name: Laboratory of Tropical Forest Ecology, Xishuangbanna Tropical Botanical Garden, Chinese Academy of Sciences, Mengla, Yunnan, China
    index: 1
  - name: Yunnan Key Laboratory of Forest Ecosystem Stability and Global Change, Xishuangbanna Tropical Botanical Garden, Chinese Academy of Sciences, Mengla, Yunnan, China
    index: 2
date: 20 December 2025
bibliography: paper.bib
---

# Summary

`calba` (Calculate Basal Area for Trees) is an R package for computing neighborhood basal area metrics for individual trees in mapped forest plots.
It calculates conspecific and total basal area within user-defined radii, supports unweighted and distance-decay kernels, and flags trees affected by plot boundaries using a simple edge correction.
The resulting outputs are model-ready covariates for neighborhood competition, growth, and survival analyses in forest ecology.

# Statement of need

Quantifying neighborhood structure is fundamental in spatial ecology and forest dynamics because competition and facilitation operate at local scales defined by distance and species identity [@Song2024a; @Hulsmann2024; @Wiegand2025].
Researchers often need to evaluate many alternative hypotheses, such as multiple radii, different decay functions, or focal subsets of species, before fitting statistical models.
Generic spatial packages can compute pairwise distances or intensity functions, but they rarely provide conspecific and total basal-area metrics with built-in edge correction and customizable kernels.
Packages such as `siplab` [@Garcia2014], `forestecology` [@Kim2021], and the ForestGEO package `fgeo` include neighborhood or spatial utilities that support a wide range of tasks.
However, their neighborhood functions are typically designed for broader plot-level workflows and therefore do not focus on producing conspecific and total basal-area summaries with flexible radii, distance-decay kernels, and simple edge treatment in a single, model-ready step.
In contrast, `calba` provides a dedicated high-performance implementation that reuses distance calculations, minimizes memory allocation, and processes radii and decay parameters in a single C++ pass.
These optimizations make `calba` substantially faster in practice, particularly when evaluating multiple radii or decay kernels, compared with naïve R-level or unoptimized C++ implementations that recompute distances repeatedly.
`calba` also returns clean, model-ready neighborhood covariates—including conspecific and total basal area, heterospecific components, decay-weighted metrics, and neighbor counts—so that researchers can proceed directly to competition, growth, or survival models without additional data wrangling.
Together, these features allow ecologists to compute comprehensive neighborhood summaries efficiently, even for large forest plots with tens or hundreds of thousands of stems.

# Usage

`neigh_ba()` is the primary interface for computing conspecific and total basal area, neighbor counts, and decay-based summaries in a single call.
The function returns a tidy data frame that can be used directly in regression and mixed-effects models.
All examples complete in under a second on a typical laptop.

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
  r_values = c(2, 5, 10)
)
```

Advanced users can access the underlying kernels through `ba_simple()`, `ba_decay()`, `count_con()`, and `count_total()`, which provide low-level components for custom workflows.

# Implementation and quality control

All core routines are implemented in `src/calba.cpp`, where distance calculations, species comparisons, and kernel applications are handled in compiled C++ loops.
The R interface provides user-friendly wrappers (`ba_simple()`, `ba_decay()`, `count_con()`, `count_total()`, `neigh_multi_r()`) with argument validation, consistent units, and shared edge-handling logic.
Automated tests using `testthat` compare all exported functions against brute-force R implementations across radii, kernels, and grouping conditions, ensuring correctness and reproducibility.
The package passes R CMD check without errors or warnings and is fully documented using roxygen2.

# Performance and workflow integration

By reusing the same pairwise distance structure for a set of radii or decay parameters, `calba` avoids repeated distance matrix recomputation that would otherwise dominate runtime when exploring model specifications.
This efficiency makes it practical to summarize neighborhood basal area for hundreds or thousands of mapped trees with dozens of kernel/radius combinations.
The outputs can be directly merged into data frames for modeling with GLMMs, GAMs, or Bayesian workflows, and optional helpers simplify counting neighbors or re-running summaries over multiple radii (`neigh_multi_r()`).

# Availability

`calba` is available on [CRAN](https://cran.r-project.org/web/packages/calba/index.html) and developed openly on [GitHub](https://github.com/mattocci27/calba).


# References
