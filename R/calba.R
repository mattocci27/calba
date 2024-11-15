#' Calculate Simple Basal Area
#'
#' This function calculates the basal area for a given set of tree diameters at breast height (DBH),
#' the number of focal trees, and a radius parameter.
#'
#' The basal area is a measure commonly used in forestry and ecology to represent the cross-sectional
#' area of a tree trunk at breast height, typically at 1.3 meters above the ground. The function
#' uses a simple formula to calculate basal area for multiple trees.
#'
#' @param dbh A numeric vector of tree diameters at breast height (DBH) in centimeters.
#' @param n_focal An integer specifying the number of focal trees to include in the calculation.
#'   This parameter limits the number of values from `dbh` that are used.
#' @param r A numeric scalar for the radius parameter. This parameter could represent a distance
#'   or scaling factor used in calculating the basal area (interpretation depends on context).
#'
#' @return A numeric vector representing the basal area for each specified tree in square centimeters.
#'   The length of the output vector is equal to `n_focal`.
#'
#' @details The formula used to calculate the basal area is based on the approximation:
#'   \deqn{\text{basal area} = \pi \cdot \left(\frac{\text{dbh}}{2}\right)^2}
#'   where DBH is the diameter at breast height.
#'
#' @examples
#' # Calculate basal area for a set of DBH values
#' sample_data <- data.frame(latin = sample(letters[1:4], 100, replace = TRUE), gx = runif(100), gy = runif(100), ba = runif(100, 10, 30))
#' ba_simple( sample_data, n_focal = 10, r = 1.5)
#'
#' @useDynLib calba, .registration = TRUE
#' @import Rcpp
#' @export
ba_simple <- function(dbh, n_focal, r) {
  calculate_basal_area_simple(dbh, n_focal, r)
}

#' Calculate Basal Area with Decay Function
#'
#' This function calculates the basal area across a given set of trees, adjusting for decay based
#' on mu values (decay parameters). It uses a zero-truncated compound Poisson-lognormal distribution
#' mixture model to handle multiple focal trees with specified decay.
#'
#' @param mu_values A numeric vector of decay parameters. Each value in `mu_values` represents
#'   a decay factor used in calculating basal area for different focal trees.
#' @param data A data frame containing the tree data, including columns such as `gx`, `gy`, and `ba`
#'   for the x-coordinates, y-coordinates, and basal area of each tree.
#' @param n_focal An integer specifying the number of focal trees to include in the calculation.
#' @param r A numeric scalar representing a decay radius or scaling factor, used to modify
#'   the decay effect in the calculations.
#'
#' @return A list with two matrices:
#' \describe{
#'   \item{`con_ba_matrix`}{A numeric matrix of basal areas with decay applied for conspecific (same species) trees.}
#'   \item{`het_ba_matrix`}{A numeric matrix of basal areas with decay applied for heterospecific (different species) trees.}
#' }
#'
#' @details This function applies a decay model to adjust basal area values based on distance
#'   and species type. The compound Poisson-lognormal distribution mixture is used to handle
#'   zero-truncation, ensuring no zero basal area values.
#'
#'   This is useful for ecological studies where the effect of tree proximity (distance-decay) and
#'   tree type (same vs. different species) are important for understanding basal area distribution.
#'
#' @examples
#' # Generate a sample data frame
#' sample_data <- data.frame(latin = sample(letters[1:4], 100, replace = TRUE), gx = runif(100), gy = runif(100), ba = runif(100, 10, 30))
#' mu_values <- c(1, 3, 5, 7)
#' ba_decay(mu_values, sample_data, n_focal = 10, r = 1.5)
#'
#' @export
ba_decay <- function(mu_values, data, n_focal, r) {
  # Calls the C++ function (assuming it’s available in the package)
  calculate_basal_area_decay(mu_values, data, n_focal, r)
}
