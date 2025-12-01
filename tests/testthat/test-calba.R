library(testthat)

brute_simple <- function(sp, gx, gy, ba, r, dist_weighted = FALSE) {
  n <- length(gx)
  con <- numeric(n)
  total <- numeric(n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      dx <- gx[j] - gx[i]
      dy <- gy[j] - gy[i]
      dist <- sqrt(dx * dx + dy * dy)
      if (dist <= 0 || dist > r) next
      contribution <- if (dist_weighted) ba[j] / dist else ba[j]
      total[i] <- total[i] + contribution
      if (sp[i] == sp[j]) {
        con[i] <- con[i] + contribution
      }
    }
  }

  list(con_ba = con, total_ba = total)
}

brute_counts <- function(sp, gx, gy, r) {
  n <- length(gx)
  con <- integer(n)
  total <- integer(n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      dx <- gx[j] - gx[i]
      dy <- gy[j] - gy[i]
      dist <- sqrt(dx * dx + dy * dy)
      if (dist > r) next
      total[i] <- total[i] + 1
      if (sp[i] == sp[j]) {
        con[i] <- con[i] + 1
      }
    }
  }

  list(con_count = con, total_count = total)
}

brute_decay <- function(mu_values, sp, gx, gy, ba, r, decay_type) {
  n <- length(gx)
  m <- length(mu_values)
  con_mat <- matrix(0, nrow = n, ncol = m)
  total_mat <- matrix(0, nrow = n, ncol = m)

  kernel <- switch(
    decay_type,
    exponential = function(ba_val, dist, mu) ba_val * exp(-dist / mu),
    `exponential-normal` = function(ba_val, dist, mu) ba_val * exp(-dist * dist / (mu * mu)),
    stop("unsupported decay type", call. = FALSE)
  )

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      dx <- gx[j] - gx[i]
      dy <- gy[j] - gy[i]
      dist <- sqrt(dx * dx + dy * dy)
      if (dist <= 0 || dist > r) next
      for (k in seq_len(m)) {
        mu <- mu_values[k]
        contribution <- kernel(ba[j], dist, mu)
        total_mat[i, k] <- total_mat[i, k] + contribution
        if (sp[i] == sp[j]) {
          con_mat[i, k] <- con_mat[i, k] + contribution
        }
      }
    }
  }

  list(con = con_mat, total = total_mat)
}

set.seed(42)
sample_data <- data.frame(
  latin = sample(letters[1:4], 100, replace = TRUE),
  gx = runif(100, 0, 10),
  gy = runif(100, 0, 10),
  ba = runif(100, 10, 30)
)
sample_sp <- sample_data$latin
sample_coords <- list(
  gx = sample_data$gx,
  gy = sample_data$gy
)
sample_ba <- sample_data$ba

test_that("ba_simple matches brute force (weighted/unweighted)", {
  res_plain <- ba_simple(sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3)
  brute_plain <- brute_simple(sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3)
  expect_equal(res_plain$con_ba, brute_plain$con_ba)
  expect_equal(res_plain$total_ba, brute_plain$total_ba)

  res_weighted <- ba_simple(sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3, dist_weighted = TRUE)
  brute_weighted <- brute_simple(sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3, dist_weighted = TRUE)
  expect_equal(res_weighted$con_ba, brute_weighted$con_ba)
  expect_equal(res_weighted$total_ba, brute_weighted$total_ba)
})

test_that("ba_decay matches brute force decay kernels", {
  mu_values <- c(1, 2)
  res_exp <- ba_decay(mu_values, sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3)
  brute_exp <- brute_decay(mu_values, sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3, "exponential")
  expect_equal(res_exp$con_ba_matrix, brute_exp$con)
  expect_equal(res_exp$total_ba_matrix, brute_exp$total)

  res_gauss <- ba_decay(mu_values, sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3, exponential_normal = TRUE)
  brute_gauss <- brute_decay(mu_values, sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3, "exponential-normal")
  expect_equal(res_gauss$con_ba_matrix, brute_gauss$con)
  expect_equal(res_gauss$total_ba_matrix, brute_gauss$total)
})

test_that("count helpers match brute force counts", {
  res_counts <- list(
    con = count_con(sample_sp, sample_coords$gx, sample_coords$gy, r = 3),
    total = count_total(sample_coords$gx, sample_coords$gy, r = 3)
  )
  brute_counts_res <- brute_counts(sample_sp, sample_coords$gx, sample_coords$gy, r = 3)
  expect_equal(res_counts$con, brute_counts_res$con_count)
  expect_equal(res_counts$total, brute_counts_res$total_count)
})

test_that("edge radii behave as expected", {
  # Tiny radius: no neighbors
  tiny <- ba_simple(sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 1e-6)
  expect_true(all(tiny$con_ba == 0))
  expect_true(all(tiny$total_ba == 0))

  counts_large <- count_total(sample_coords$gx, sample_coords$gy, r = 100)
  expect_equal(counts_large, rep(length(sample_sp) - 1L, length(sample_sp)))
})

test_that("edge_correction 'safe' skips edge focal trees", {
  r_val <- 1
  safe_res <- ba_simple(
    sample_sp,
    sample_coords$gx,
    sample_coords$gy,
    sample_ba,
    r = r_val,
    edge_correction = "safe"
  )
  expect_true(any(is.na(safe_res$con_ba)))
  expect_true(any(is.na(safe_res$total_ba)))

  safe_counts <- count_total(
    sample_coords$gx,
    sample_coords$gy,
    r = r_val,
    edge_correction = "safe"
  )
  expect_true(any(is.na(safe_counts)))

  decay_safe <- ba_decay(
    c(1, 2),
    sample_sp,
    sample_coords$gx,
    sample_coords$gy,
    sample_ba,
    r = r_val,
    edge_correction = "safe"
  )
  expect_true(any(is.na(decay_safe$con_ba_matrix)))
})

test_that("user-supplied bounds control safe edge detection", {
  r_val <- 1
  xmin <- min(sample_coords$gx)
  xmax <- max(sample_coords$gx)
  ymin <- min(sample_coords$gy)
  ymax <- max(sample_coords$gy)
  expected_edge <- !(
    (sample_coords$gx - r_val >= xmin) &
      (sample_coords$gx + r_val <= xmax) &
      (sample_coords$gy - r_val >= ymin) &
      (sample_coords$gy + r_val <= ymax)
  )

  default_safe <- ba_simple(
    sample_sp,
    sample_coords$gx,
    sample_coords$gy,
    sample_ba,
    r = r_val,
    edge_correction = "safe"
  )
  expect_equal(is.na(default_safe$con_ba), expected_edge)
  custom_bounds <- c(-1, 11, -1, 11)
  bounded_safe <- ba_simple(
    sample_sp,
    sample_coords$gx,
    sample_coords$gy,
    sample_ba,
    r = r_val,
    edge_correction = "safe",
    bounds = custom_bounds
  )
  expect_false(any(is.na(bounded_safe$con_ba)))
  expect_true(all(sample_coords$gx <= custom_bounds[2]))
})

test_that("neigh_multi_r matches repeated ba_simple", {
  radii <- c(1, 3)
  res_multi <- neigh_multi_r(sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, radii)
  expect_equal(res_multi$radius, rep(radii, each = length(sample_sp)))
  expect_equal(
    res_multi$con_ba[res_multi$radius == radii[2]],
    ba_simple(sample_sp, sample_coords$gx, sample_coords$gy, sample_ba, r = 3)$con_ba
  )
})
