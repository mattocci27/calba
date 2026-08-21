library(testthat)

brute_simple <- function(sp, gx, gy, ba, r, dist_weighted = FALSE,
                         zero_distance = c("include", "exclude")) {
  zero_distance <- match.arg(zero_distance)
  n <- length(gx)
  con <- numeric(n)
  total <- numeric(n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      dx <- gx[j] - gx[i]
      dy <- gy[j] - gy[i]
      dist <- sqrt(dx * dx + dy * dy)
      if (dist > r || (dist == 0 && zero_distance == "exclude")) next
      if (dist == 0 && dist_weighted) {
        stop("zero distance is undefined for inverse-distance weighting", call. = FALSE)
      }
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

brute_decay <- function(mu_values, sp, gx, gy, ba, r, decay_type,
                        zero_distance = c("include", "exclude")) {
  zero_distance <- match.arg(zero_distance)
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
      if (dist > r || (dist == 0 && zero_distance == "exclude")) next
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

test_that("zero-distance policy is explicit and exponential kernels can include co-located stems", {
  sp <- c("a", "a", "b")
  gx <- c(0, 0, 1)
  gy <- c(0, 0, 0)
  ba <- c(1, 2, 3)
  r <- 2

  expect_error(
    ba_decay(1, sp, gx, gy, ba, r),
    "identical coordinates"
  )
  expect_error(
    count_total(gx, gy, r),
    "identical coordinates"
  )

  expected_simple <- list(con_ba = c(2, 1, 0), total_ba = c(5, 4, 3))
  simple_include <- ba_simple(
    sp, gx, gy, ba, r,
    zero_distance = "include"
  )
  expect_equal(simple_include, expected_simple)
  expect_equal(
    count_con(sp, gx, gy, r, zero_distance = "include"),
    c(1, 1, 0)
  )
  expect_equal(
    count_total(gx, gy, r, zero_distance = "include"),
    c(2, 2, 2)
  )

  exp_weight <- exp(-1 / 2)
  expn_weight <- exp(-1 / 4)
  expected_exp_total <- matrix(c(2 + 3 * exp_weight, 1 + 3 * exp_weight,
                                  3 * exp_weight), ncol = 1)
  expected_expn_total <- matrix(c(2 + 3 * expn_weight, 1 + 3 * expn_weight,
                                   3 * expn_weight), ncol = 1)
  expected_con <- matrix(c(2, 1, 0), ncol = 1)

  exp_include <- ba_decay(2, sp, gx, gy, ba, r, zero_distance = "include")
  expect_equal(exp_include$con_ba_matrix, expected_con)
  expect_equal(exp_include$total_ba_matrix, expected_exp_total)

  expn_include <- ba_decay(
    2, sp, gx, gy, ba, r,
    exponential_normal = TRUE,
    zero_distance = "include"
  )
  expect_equal(expn_include$con_ba_matrix, expected_con)
  expect_equal(expn_include$total_ba_matrix, expected_expn_total)

  neighborhood_include <- neigh_ba(
    sp, gx, gy, ba, r,
    mu_values = 2,
    zero_distance = "include"
  )
  expect_equal(neighborhood_include$summary$total_ba, expected_simple$total_ba)
  expect_equal(neighborhood_include$decay$total_ba, as.vector(expected_exp_total))

  multi_include <- neigh_multi_r(
    sp, gx, gy, ba, r_values = r,
    zero_distance = "include"
  )
  expect_equal(multi_include$total_ba, expected_simple$total_ba)

  simple_exclude <- ba_simple(
    sp, gx, gy, ba, r,
    zero_distance = "exclude"
  )
  expect_equal(simple_exclude$total_ba, c(3, 3, 3))
  expect_error(
    ba_simple(sp, gx, gy, ba, r, dist_weighted = TRUE, zero_distance = "include"),
    "undefined for inverse-distance weighting"
  )
})

test_that("duplicate_coordinates reports but does not modify coordinate records", {
  duplicates <- duplicate_coordinates(c(0, 0, 1, 2, 2), c(0, 0, 1, 2, 2))
  expect_equal(
    duplicates,
    data.frame(
      tree_id = c(1L, 2L, 4L, 5L),
      gx = c(0, 0, 2, 2),
      gy = c(0, 0, 2, 2),
      coordinate_group = c(1L, 1L, 2L, 2L)
    )
  )
  expect_equal(
    duplicate_coordinates(c(0, 1), c(0, 1)),
    data.frame(
      tree_id = integer(), gx = numeric(), gy = numeric(), coordinate_group = integer()
    )
  )
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
