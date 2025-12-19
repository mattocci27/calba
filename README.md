# calba: Efficient neighborhood basal area metrics for trees  <img src="man/figures/icon.png" width = "200" align="right" />

`calba` stands for *Calculate Basal Area for Trees*.
It provides fast, Rcpp-backed routines to generate conspecific and total neighborhood basal area metrics for every stem in a mapped forest plot.
The package keeps radius, decay kernel, and edge-correction choices explicit so researchers can explore alternative neighborhood hypotheses without re-implementing the costly distance calculations in R.

## Installation

Install the development version directly from GitHub:

```r
remotes::install_github("mattocci27/calba")
```

## Basic usage

Below is a reproducible workflow that mirrors a typical research use case:
compute neighborhood basal area under multiple radii and kernels, then examine how those covariates might enter a model.

```r
library(calba)

set.seed(42)
sample_trees <- data.frame(
  sp = sample(c("oak", "pine", "birch"), 40, replace = TRUE),
  gx = runif(40, 0, 20),
  gy = runif(40, 0, 20),
  ba = runif(40, 5, 35)
)

# conspecific/total basal area within r = 5 units
neighborhood <- ba_simple(
  sp = sample_trees$sp,
  gx = sample_trees$gx,
  gy = sample_trees$gy,
  ba = sample_trees$ba,
  r = 5,
  dist_weighted = TRUE,
  edge_correction = "safe"
)

head(neighborhood)
```

Call `ba_decay()` when you want to compare exponential and exponential-normal distance decay kernels or evaluate a grid of decay parameters (`mu_values`). Use `count_con()`, `count_total()`, and `neigh_multi_r()` to extract counts or re-run the same calculations over multiple radii without recomputing the distance matrix.

## Example research scenario

Researchers studying competition can use `calba` to generate covariates for growth or survival models. A typical pipeline looks like:

1. Map tree locations/DBH and compute basal area.
2. Use `ba_simple()` and `ba_decay()` to summarize neighborhood structure for a set of biologically motivated radii and kernels.
3. Feed the resulting matrices (or selected radius/kern) into GLMMs or Bayesian models to test hypotheses about conspecific vs. heterospecific competition.

Because the heavy lifting is done in C++, the same distances are reused when evaluating multiple radii or kernels, which keeps workflows feasible on large datasets.

## Testing and checks

Run the package tests with:

```r
devtools::test()
```

and validate the package before submission with:

```r
devtools::check()
```
