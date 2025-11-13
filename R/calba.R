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
ba_simple <- function(sp, gx, gy, ba, r, dist_weighted = FALSE) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))

  calculate_basal_area_simple(sp, gx, gy, ba, r, dist_weighted)
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
ba_decay <- function(mu_values, sp, gx, gy, ba, r, exponential_normal = FALSE) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))
  mu_values <- validate_mu_values(mu_values)

  decay_type <- if (exponential_normal) "exponential-normal" else "exponential"
  calculate_basal_area_decay(mu_values, sp, gx, gy, ba, r, decay_type)
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
count_con <- function(sp, gx, gy, r) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))

  count_con_cpp(sp, gx, gy, r)
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
count_total <- function(gx, gy, r) {
  validate_xy(gx, gy, r)

  count_total_cpp(gx, gy, r)
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
#'
#' @return A list with
#' * `summary`: data frame with `tree_id`, `species`, `con_ba`, `total_ba`, `con_count`, `total_count`.
#' * `decay`: (`NULL` or) long data frame with `tree_id`, `species`, `mu`, `con_ba`, `total_ba`.
#'
#' @examples
#' sample_data <- data.frame(
#'   latin = sample(letters[1:4], 10, replace = TRUE),
#'   gx = runif(10, 0, 10),
#'   gy = runif(10, 0, 10),
#'   ba = runif(10, 10, 30)
#' )
#' neighborhood_ba(
#'   sp = sample_data$latin,
#'   gx = sample_data$gx,
#'   gy = sample_data$gy,
#'   ba = sample_data$ba,
#'   r = 3,
#'   mu_values = c(1, 3)
#' )
#' @export
neighborhood_ba <- function(sp, gx, gy, ba, r,
                            mu_values = NULL,
                            dist_weighted = FALSE,
                            exponential_normal = FALSE) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))
  if (!is.null(mu_values)) {
    mu_values <- validate_mu_values(mu_values)
  }

  ba_out <- ba_simple(sp, gx, gy, ba, r, dist_weighted = dist_weighted)
  summary_tbl <- data.frame(
    tree_id = seq_len(length(gx)),
    species = sp,
    con_ba = ba_out$con_ba,
    total_ba = ba_out$total_ba,
    con_count = count_con(sp, gx, gy, r),
    total_count = count_total(gx, gy, r),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  decay_tbl <- NULL
  if (!is.null(mu_values)) {
    decay_tbl <- ba_decay_long(
      mu_values = mu_values,
      sp = sp,
      gx = gx,
      gy = gy,
      ba = ba,
      r = r,
      exponential_normal = exponential_normal
    )
  }

  list(
    summary = summary_tbl,
    decay = decay_tbl
  )
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
#'
#' @return A data frame with `tree_id`, `species`, `mu`, `con_ba`, and `total_ba`.
#'
#' @export
ba_decay_long <- function(mu_values, sp, gx, gy, ba, r, exponential_normal = FALSE) {
  validate_xy(gx, gy, r)
  sp <- validate_species(sp, length(gx))
  ba <- validate_ba(ba, length(gx))
  mu_values <- validate_mu_values(mu_values)

  decay_res <- ba_decay(
    mu_values = mu_values,
    sp = sp,
    gx = gx,
    gy = gy,
    ba = ba,
    r = r,
    exponential_normal = exponential_normal
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

validate_xy <- function(gx, gy, r) {
  if (!is.numeric(gx) || !is.numeric(gy)) {
    stop("`gx` and `gy` must be numeric vectors", call. = FALSE)
  }
  if (length(gx) != length(gy)) {
    stop("`gx` and `gy` must have the same length", call. = FALSE)
  }
  if (length(r) != 1 || !is.finite(r) || r <= 0) {
    stop("`r` must be a single finite, positive number", call. = FALSE)
  }
  if (any(!is.finite(gx)) || any(!is.finite(gy))) {
    stop("`gx` and `gy` must contain only finite values", call. = FALSE)
  }
  invisible(NULL)
}

validate_species <- function(sp, n) {
  sp <- as.character(sp)
  if (length(sp) != n) {
    stop("`sp` must have the same length as `gx`/`gy`", call. = FALSE)
  }
  if (any(is.na(sp))) {
    stop("`sp` cannot contain `NA` values", call. = FALSE)
  }
  sp
}

validate_ba <- function(ba, n) {
  if (!is.numeric(ba)) {
    stop("`ba` must be numeric", call. = FALSE)
  }
  if (length(ba) != n) {
    stop("`ba` must have the same length as the coordinate vectors", call. = FALSE)
  }
  if (any(!is.finite(ba))) {
    stop("`ba` must contain only finite values", call. = FALSE)
  }
  ba
}

validate_mu_values <- function(mu_values) {
  if (!is.numeric(mu_values)) {
    stop("`mu_values` must be numeric", call. = FALSE)
  }
  if (length(mu_values) == 0) {
    stop("`mu_values` must contain at least one value", call. = FALSE)
  }
  if (any(!is.finite(mu_values)) || any(mu_values <= 0)) {
    stop("`mu_values` must be positive and finite", call. = FALSE)
  }
  mu_values
}
