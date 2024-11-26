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
#'   r = 3
#' )
#'
#' @export
ba_decay <- function(mu_values, sp, gx, gy, ba, r) {
  calculate_basal_area_decay(mu_values, sp, gx, gy, ba, r)
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
  count_total_cpp(gx, gy, r)
}
