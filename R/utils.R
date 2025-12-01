## Internal input-validation helpers used across the exported R wrappers.
## Keep these light and fail fast so that the R layer never sends bad data to the C++ core.

## Ensure coordinates and radius are numeric, finite, and length-compatible.
## This protects each exported wrapper before it reaches the C++ layer.
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

## Normalize species labels and reject mismatched lengths/NA entries.
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

## Basal-area vector must align with tree coordinates and be numeric/finite.
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

## Decay parameters (`mu`) must be strictly positive for kernel smoothing.
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

## Multi-radius helpers rely on a finite positive radius vector.
validate_r_values <- function(r_values) {
  if (!is.numeric(r_values)) {
    stop("`r_values` must be numeric", call. = FALSE)
  }
  if (length(r_values) == 0) {
    stop("`r_values` must contain at least one value", call. = FALSE)
  }
  if (any(!is.finite(r_values)) || any(r_values <= 0)) {
    stop("`r_values` must be positive and finite", call. = FALSE)
  }
  r_values
}

## Optional user-supplied plot bounds: NULL or c(xmin, xmax, ymin, ymax).
## Ensures inputs are finite, ordered, and cover all coordinates.
validate_bounds <- function(bounds, gx, gy) {
  if (is.null(bounds)) {
    return(NULL)
  }
  if (!is.numeric(bounds) || length(bounds) != 4) {
    stop("`bounds` must be a numeric vector of length 4: xmin, xmax, ymin, ymax", call. = FALSE)
  }
  if (any(!is.finite(bounds))) {
    stop("`bounds` must be finite", call. = FALSE)
  }
  xmin <- bounds[[1]]
  xmax <- bounds[[2]]
  ymin <- bounds[[3]]
  ymax <- bounds[[4]]
  if (!(xmin < xmax && ymin < ymax)) {
    stop("`bounds` must satisfy xmin < xmax and ymin < ymax", call. = FALSE)
  }
  if (any(gx < xmin | gx > xmax | gy < ymin | gy > ymax)) {
    stop("All coordinates must fall within `bounds`", call. = FALSE)
  }
  bounds
}
