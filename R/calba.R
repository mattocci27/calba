#' Calculate Simple Basal Area
#'
#' This function calculates the total basal area (conspecific + heterospecific) for a given set of tree data,
#' focusing on the number of focal trees and a radius parameter. The basal area represents the cumulative 
#' cross-sectional area of tree trunks adjusted for a simple decay model.
#'
#' @param data A data frame containing the tree data, including columns such as `latin` (species name), 
#'   `gx`, `gy` (coordinates), and `ba` (basal area).
#' @param n_focal An integer specifying the number of focal trees to include in the calculation.
#' @param r A numeric scalar for the radius parameter. This parameter represents the maximum distance 
#'   to consider when summing basal areas of neighboring trees.
#'
#' @return A list with two numeric vectors:
#' \describe{
#'   \item{`con_ba`}{A numeric vector representing the cumulative basal area of conspecific trees within the radius.}
#'   \item{`total_ba`}{A numeric vector representing the cumulative basal area of all trees (conspecific + heterospecific) within the radius.}
#' }
#'
#' @details The function applies a simple decay model where the basal area contribution of each tree decreases
#'   proportionally with its distance from the focal tree:
#'   \deqn{\text{decayed basal area} = \frac{\text{ba}}{\text{dist}}}
#'   where `dist` is the Euclidean distance from the focal tree to its neighbors.
#'
#' @examples
#' # Generate a sample dataset
#' sample_data <- data.frame(
#'   latin = sample(letters[1:4], 100, replace = TRUE),
#'   gx = runif(100),
#'   gy = runif(100),
#'   ba = runif(100, 10, 30)
#' )
#'
#' ba_simple(sample_data, n_focal = 10, r = 1.5)
#'
#' @useDynLib calba, .registration = TRUE
#' @import Rcpp
#' @export
ba_simple <- function(data, n_focal, r) {
  calculate_basal_area_simple(data, n_focal, r)
}

#' Calculate Basal Area with Decay Function
#'
#' This function calculates the basal area across a given set of trees, applying a decay effect
#' based on distance and species identity for each tree within a given radius.
#'
#' @param mu_values A numeric vector of decay parameters. Each value in `mu_values` represents
#'   a decay factor that modifies how the basal area contribution diminishes with distance.
#' @param data A data frame containing the tree data, including columns such as `latin` (species name), 
#'   `gx`, `gy` (coordinates), and `ba` (basal area).
#' @param n_focal An integer specifying the number of focal trees to include in the calculation.
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
#'   This approach allows the analysis of how tree proximity (distance-decay) and species identity 
#'   affect the distribution of basal area in an ecological setting.
#'
#' @examples
#' # Generate a sample dataset
#' sample_data <- data.frame(
#'   latin = sample(letters[1:4], 100, replace = TRUE),
#'   gx = runif(100),
#'   gy = runif(100),
#'   ba = runif(100, 10, 30)
#' )
#' mu_values <- c(1, 3, 5, 7)
#' ba_decay(mu_values, sample_data, n_focal = 10, r = 1.5)
#'
#' @export
ba_decay <- function(mu_values, data, n_focal, r) {
  calculate_basal_area_decay(mu_values, data, n_focal, r)
}


