#' Calculate Simple Basal Area
#'
#' This function calculates the total basal area (conspecific + heterospecific) for a given set of tree data,
#' focusing on the number of focal trees and a radius parameter. The basal area represents the cumulative
#' cross-sectional area of tree trunks, either unadjusted or adjusted for a simple distance-weighted decay model.
#'
#' @param sp A character vector containing species names for each tree.
#' @param gx A numeric vector of x-coordinates for the trees (e.g., 0 to 10 for realistic spatial data).
#' @param gy A numeric vector of y-coordinates for the trees (e.g., 0 to 10 for realistic spatial data).
#' @param ba A numeric vector of basal area values for the trees (e.g., cross-sectional area at breast height).
#' @param r A numeric scalar for the radius parameter. This parameter represents the maximum distance
#'   to consider when summing basal areas of neighboring trees (e.g., 1–5 units for spatially distributed data).
#' @param dist_weighted Logical. If `TRUE`, use the distance-weighted approach (`ba / dist`);
#'   if `FALSE`, use the unadjusted `ba`.
#' @param edge_correction Character. Use `"none"` (default) to return all focal trees or `"safe"`
#'   to skip neighborhood calculations for trees within `r` of the plot edges and mark their output
#'   as `NA`.
#' @param bounds Optional numeric vector `c(xmin, xmax, ymin, ymax)` giving the plot extent.
#'   When `NULL`, the range of `gx`/`gy` is used; supply bounds if your data do not span the full plot.
#' 
#' @return A list with two numeric vectors:
#' \describe{
#'   \item{`con_ba`}{A numeric vector representing the cumulative basal area of conspecific trees within the radius.}
#'   \item{`total_ba`}{A numeric vector representing the cumulative basal area of all trees (conspecific + heterospecific) within the radius.}
#' }
#'
#' @details The function either:
#'   - Calculates the unadjusted basal area if `dist_weighted = FALSE`.
#'   - Applies a simple distance-weighted approach (`ba / dist`) if `dist_weighted = TRUE`.
#'
#' @examples
#' # Generate a sample dataset
#' set.seed(42)  # For reproducibility
#' sample_data <- data.frame(
#'   latin = sample(letters[1:4], 100, replace = TRUE),
#'   gx = runif(100, 0, 10),  # Spatial coordinates between 0 and 10
#'   gy = runif(100, 0, 10),
#'   ba = runif(100, 10, 30)  # Basal area between 10 and 30
#' )
#'
#' # Calculate with distance weighting
#' ba_simple(
#'   sp = sample_data$latin,
#'   gx = sample_data$gx,
#'   gy = sample_data$gy,
#'   ba = sample_data$ba,
#'   r = 3,  # Radius within the spatial scale
#'   dist_weighted = TRUE
#' )
#'
#' # Calculate without distance weighting
#' ba_simple(
#'   sp = sample_data$latin,
#'   gx = sample_data$gx,
#'   gy = sample_data$gy,
#'   ba = sample_data$ba,
#'   r = 3,
#'   dist_weighted = FALSE
#' )
#'
#' @useDynLib calba, .registration = TRUE
#' @import Rcpp
#' @export
ba_simple <- function(sp, gx, gy, ba, r, dist_weighted = FALSE,
                      edge_correction = c("none", "safe"),
                      bounds = NULL) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))
  bounds <- validate_bounds(bounds, gx, gy)

  edge_correction <- match.arg(edge_correction)

  calculate_basal_area_simple(
    sp, gx, gy, ba, r,
    dist_weighted = dist_weighted,
    edge_correction = edge_correction,
    bounds = bounds
  )
}

#' Calculate Basal Area with Decay Function
#'
#' This function calculates the basal area across a given set of trees, applying a decay effect
#' based on distance and species identity for each tree within a given radius.
#'
#' @param mu_values A numeric vector of decay parameters. Each value in `mu_values` represents
#'   a decay factor that modifies how the basal area contribution diminishes with distance.
#' @param sp A character vector containing species names for each tree.
#' @param gx A numeric vector of x-coordinates for the trees.
#' @param gy A numeric vector of y-coordinates for the trees.
#' @param ba A numeric vector of basal area values for the trees.
#' @param r A numeric scalar representing the radius to consider for neighboring trees.
#' @param exponential_normal A logical value. If `FALSE` (default), use exponential decay.
#' If `TRUE`, use exponential-normal decay.
#' @param edge_correction Character. See `ba_simple()` for the `"safe"` behavior that skips edge trees.
#' @param bounds Optional numeric vector `c(xmin, xmax, ymin, ymax)` giving the plot extent.
#'   When `NULL`, the range of `gx`/`gy` is used; supply bounds if your data do not span the full plot.
#'
#' @return A list with two matrices:
#' \describe{
#'   \item{`con_ba_matrix`}{A numeric matrix of basal areas with decay applied for conspecific (same species) trees.}
#'   \item{`total_ba_matrix`}{A numeric matrix of basal areas with decay applied for all trees (conspecific + heterospecific).}
#' }
#'
#' @details The function applies an exponential decay model where the basal area contribution
#'   diminishes with distance from the focal tree:
#'   \deqn{\text{decayed basal area} = \text{ba} \cdot \exp\left(-\frac{\text{dist}}{\mu}\right)}
#'   where `mu` is the decay parameter, `ba` is the basal area, and `dist` is the Euclidean distance
#'   between trees.
#'
#' @examples
#' # Generate a sample dataset
#' set.seed(42)
#' sample_data <- data.frame(
#'   latin = sample(letters[1:4], 100, replace = TRUE),
#'   gx = runif(100, 0, 10),
#'   gy = runif(100, 0, 10),
#'   ba = runif(100, 10, 30)
#' )
#' mu_values <- c(1, 3, 5, 7)
#' ba_decay(
#'   mu_values = mu_values,
#'   sp = sample_data$latin,
#'   gx = sample_data$gx,
#'   gy = sample_data$gy,
#'   ba = sample_data$ba,
#'   r = 3,
#'   exponential_normal = FALSE
#' )
#'
#' @export
ba_decay <- function(mu_values, sp, gx, gy, ba, r, exponential_normal = FALSE,
                     edge_correction = c("none", "safe"),
                     bounds = NULL) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))
  mu_values <- validate_mu_values(mu_values)
  bounds <- validate_bounds(bounds, gx, gy)

  edge_correction <- match.arg(edge_correction)

  decay_type <- if (exponential_normal) "exponential-normal" else "exponential"
  calculate_basal_area_decay(
    mu_values, sp, gx, gy, ba, r, decay_type,
    edge_correction = edge_correction,
    bounds = bounds
  )
}

#' Count Conspecific Trees
#'
#' This function counts the number of conspecific trees within a given radius for each focal tree.
#'
#' @param sp A character vector of species names.
#' @param gx A numeric vector of x-coordinates for the trees.
#' @param gy A numeric vector of y-coordinates for the trees.
#' @param r A numeric scalar for the radius parameter.
#' 
#' @return A numeric vector containing the count of conspecific trees within the radius for each focal tree.
#' @param edge_correction Character; see `ba_simple()` for the `"safe"` option that skips focal trees near the boundary.
#' @param bounds Optional numeric vector `c(xmin, xmax, ymin, ymax)` giving the plot extent.
#'   When `NULL`, the range of `gx`/`gy` is used; supply bounds if your data do not span the full plot.
#'
#' @examples
#' sample_data <- data.frame(
#'   latin = sample(letters[1:4], 100, replace = TRUE),
#'   gx = runif(100, 0, 10),
#'   gy = runif(100, 0, 10)
#' )
#' count_con(
#'   sp = sample_data$latin,
#'   gx = sample_data$gx,
#'   gy = sample_data$gy,
#'   r = 3
#' )
#'
#' @export
count_con <- function(sp, gx, gy, r, edge_correction = c("none", "safe"), bounds = NULL) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  bounds <- validate_bounds(bounds, gx, gy)

  edge_correction <- match.arg(edge_correction)

  count_con_cpp(sp, gx, gy, r,
                edge_correction = edge_correction,
                bounds = bounds)
}

#' Count Total Trees
#'
#' This function counts the total number of trees within a given radius for each focal tree.
#'
#' @param gx A numeric vector of x-coordinates for the trees.
#' @param gy A numeric vector of y-coordinates for the trees.
#' @param r A numeric scalar for the radius parameter.
#' 
#' @return A numeric vector containing the count of all trees within the radius for each focal tree.
#' @param edge_correction Character; see `ba_simple()` for the `"safe"` option that skips focal trees close to the edges.
#' @param bounds Optional numeric vector `c(xmin, xmax, ymin, ymax)` giving the plot extent.
#'   When `NULL`, the range of `gx`/`gy` is used; supply bounds if your data do not span the full plot.
#'
#' @examples
#' sample_data <- data.frame(
#'   gx = runif(100, 0, 10),
#'   gy = runif(100, 0, 10)
#' )
#' count_total(
#'   gx = sample_data$gx,
#'   gy = sample_data$gy,
#'   r = 3
#' )
#'
#' @export
count_total <- function(gx, gy, r, edge_correction = c("none", "safe"), bounds = NULL) {
  validate_xy(gx, gy, r)
  bounds <- validate_bounds(bounds, gx, gy)

  edge_correction <- match.arg(edge_correction)
  count_total_cpp(gx, gy, r,
                  edge_correction = edge_correction,
                  bounds = bounds)
}

#' Neighborhood summaries for basal area and counts
#'
#' Provide a tidy summary of pairwise neighborhood basal area and counts,
#' optionally including decay results for a vector of decay parameters.
#'
#' @param sp A character or factor vector of species names.
#' @param gx Numeric x-coordinates of each tree.
#' @param gy Numeric y-coordinates of each tree.
#' @param ba Numeric basal area of each tree.
#' @param r Positive numeric radius within which neighbors are considered.
#' @param mu_values Optional numeric vector of decay parameters. When `NULL`, the decay table is omitted.
#' @param dist_weighted Logical flag passed to `ba_simple` to use a simple `ba / dist` weighting when `TRUE`.
#' @param exponential_normal Logical passed to `ba_decay` to select the exponential-normal kernel.
#' @param edge_correction Character; see `ba_simple()` for the `"safe"` option that skips edge focal trees.
#' @param bounds Optional numeric vector `c(xmin, xmax, ymin, ymax)` giving the plot extent.
#'   When `NULL`, the range of `gx`/`gy` is used; supply bounds if your data do not span the full plot.
#'
#' @return A list with
#' * `summary`: data frame with `tree_id`, `species`, `con_ba`, `total_ba`, `con_count`, `total_count`.
#' * `decay`: (`NULL` or) long data frame with `tree_id`, `species`, `mu`, `con_ba`, `total_ba`.
#'
#' The `summary` component also includes derived columns:
#' `prop_con_ba`, `het_ba`, `het_count`, and `competition_index`.
#' @examples
#' sample_data <- data.frame(
#'   latin = sample(letters[1:4], 10, replace = TRUE),
#'   gx = runif(10, 0, 10),
#'   gy = runif(10, 0, 10),
#'   ba = runif(10, 10, 30)
#' )
#' neigh_ba(
#'   sp = sample_data$latin,
#'   gx = sample_data$gx,
#'   gy = sample_data$gy,
#'   ba = sample_data$ba,
#'   r = 3,
#'   mu_values = c(1, 3)
#' )
#' @name neigh_ba
#' @rdname neigh_ba
#' @export
neigh_ba <- function(sp, gx, gy, ba, r,
                     mu_values = NULL,
                     dist_weighted = FALSE,
                     exponential_normal = FALSE,
                     edge_correction = c("none", "safe"),
                     bounds = NULL) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))
  bounds <- validate_bounds(bounds, gx, gy)
  edge_correction <- match.arg(edge_correction)
  if (!is.null(mu_values)) {
    mu_values <- validate_mu_values(mu_values)
  }

  ba_out <- ba_simple(
    sp, gx, gy, ba, r,
    dist_weighted = dist_weighted,
    edge_correction = edge_correction,
    bounds = bounds
  )
  summary_tbl <- data.frame(
    tree_id = seq_len(length(gx)),
    species = sp,
    con_ba = ba_out$con_ba,
    total_ba = ba_out$total_ba,
    con_count = count_con(sp, gx, gy, r, edge_correction = edge_correction, bounds = bounds),
    total_count = count_total(gx, gy, r, edge_correction = edge_correction, bounds = bounds),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  summary_tbl <- add_derived_neighborhood_metrics(summary_tbl)

  decay_tbl <- NULL
  if (!is.null(mu_values)) {
    decay_tbl <- ba_decay_long(
      mu_values = mu_values,
      sp = sp,
      gx = gx,
      gy = gy,
      ba = ba,
      r = r,
      exponential_normal = exponential_normal,
      edge_correction = edge_correction,
      bounds = bounds
    )
  }

  list(
    summary = summary_tbl,
    decay = decay_tbl
  )
}

#' Neighborhood metrics across multiple radii
#'
#' Compute basal area sums and counts for each tree across multiple radii in one pass.
#'
#' @param sp A character or factor vector of species names.
#' @param gx Numeric x-coordinates of the trees.
#' @param gy Numeric y-coordinates of the trees.
#' @param ba Numeric basal area for each tree.
#' @param r_values Numeric vector of radii to evaluate.
#' @param dist_weighted Logical flag to use `ba / dist` weighting within each radius.
#' @param edge_correction Character; see `ba_simple()` for the `"safe"` option.
#' @param bounds Optional numeric vector `c(xmin, xmax, ymin, ymax)` giving the plot extent.
#'   When `NULL`, the range of `gx`/`gy` is used; supply bounds if your data do not span the full plot.
#'
#' @return A tidy tibble with `tree_id`, `species`, `radius`, `con_ba`, `total_ba`,
#' `con_count`, `total_count`, `prop_con_ba`, `het_ba`, `het_count`, and
#' `competition_index`.
#'
#' @examples
#' sample_data <- data.frame(
#'   latin = sample(letters[1:4], 10, replace = TRUE),
#'   gx = runif(10, 0, 10),
#'   gy = runif(10, 0, 10),
#'   ba = runif(10, 10, 30)
#' )
#' neigh_multi_r(
#'   sp = sample_data$latin,
#'   gx = sample_data$gx,
#'   gy = sample_data$gy,
#'   ba = sample_data$ba,
#'   r_values = c(3, 5)
#' )
#'
#' @name neigh_multi_r
#' @rdname neigh_multi_r
#' @export
neigh_multi_r <- function(sp, gx, gy, ba, r_values, dist_weighted = FALSE,
                          edge_correction = c("none", "safe"),
                          bounds = NULL) {
  r_values <- validate_r_values(r_values)
  max_r <- max(r_values)

  validate_xy(gx, gy, max_r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))
  bounds <- validate_bounds(bounds, gx, gy)

  edge_correction <- match.arg(edge_correction)
  res <- calculate_neighborhood_multi_radius(
    sp, gx, gy, ba, r_values,
    dist_weighted = dist_weighted,
    edge_correction = edge_correction,
    bounds = bounds
  )

  n <- length(gx)
  df <- data.frame(
    tree_id = rep(seq_len(n), times = length(r_values)),
    species = rep(sp, times = length(r_values)),
    radius = rep(r_values, each = n),
    con_ba = as.vector(res$con_ba),
    total_ba = as.vector(res$total_ba),
    con_count = as.vector(res$con_count),
    total_count = as.vector(res$total_count),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  df <- add_derived_neighborhood_metrics(df)

  df
}

add_derived_neighborhood_metrics <- function(summary_tbl) {
  required <- c("con_ba", "total_ba", "con_count", "total_count")
  if (!all(required %in% names(summary_tbl))) {
    stop(
      "`summary_tbl` must contain columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }

  total_ba <- summary_tbl$total_ba
  total_count <- summary_tbl$total_count

  summary_tbl$prop_con_ba <- ifelse(total_ba == 0, 0, summary_tbl$con_ba / total_ba)
  summary_tbl$het_ba <- total_ba - summary_tbl$con_ba
  summary_tbl$het_count <- total_count - summary_tbl$con_count
  summary_tbl$competition_index <- ifelse(total_count > 0, total_ba / total_count, 0)

  summary_tbl
}

#' Long-format decay table
#'
#' Transform the matrices returned by `ba_decay` into a tidy table that can be
#' joined with other data or mapped with `ggplot2`.
#'
#' @param mu_values Numeric vector of decay parameters.
#' @param sp Character/factor species vector.
#' @param gx Numeric x-coordinates.
#' @param gy Numeric y-coordinates.
#' @param ba Numeric basal area.
#' @param r Positive radius threshold.
#' @param exponential_normal Logical passed to `ba_decay`.
#' @param edge_correction Character; see `ba_simple()` for the `"safe"` option.
#' @param bounds Optional numeric vector `c(xmin, xmax, ymin, ymax)` giving the plot extent.
#'   When `NULL`, the range of `gx`/`gy` is used; supply bounds if your data do not span the full plot.
#'
#' @return A data frame with `tree_id`, `species`, `mu`, `con_ba`, and `total_ba`.
#'
#' @export
ba_decay_long <- function(mu_values, sp, gx, gy, ba, r, exponential_normal = FALSE,
                          edge_correction = c("none", "safe"),
                          bounds = NULL) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))
  mu_values <- validate_mu_values(mu_values)
  bounds <- validate_bounds(bounds, gx, gy)

  edge_correction <- match.arg(edge_correction)

  decay_res <- ba_decay(
    mu_values = mu_values,
    sp = sp,
    gx = gx,
    gy = gy,
    ba = ba,
    r = r,
    exponential_normal = exponential_normal,
    edge_correction = edge_correction,
    bounds = bounds
  )

  n <- length(gx)
  mu_rep <- rep(mu_values, each = n)
  data.frame(
    tree_id = rep(seq_len(n), times = length(mu_values)),
    species = rep(sp, times = length(mu_values)),
    mu = mu_rep,
    con_ba = as.vector(decay_res$con_ba_matrix),
    total_ba = as.vector(decay_res$total_ba_matrix),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
