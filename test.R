
set.seed(42)

sample_data <- data.frame(
   latin = sample(letters[1:4], 10, replace = TRUE),
   gx = runif(10, 0, 10),
   gy = runif(10, 0, 10),
   ba = runif(10, 10, 30)
 )

neigh_ba(
  sp = sample_data$latin,
  gx = sample_data$gx,
  gy = sample_data$gy,
  ba = sample_data$ba,
  r = 3,
  mu_values = c(1, 3)
)

ba_simple(
  sp = sample_data$latin,
  gx = sample_data$gx,
  gy = sample_data$gy,
  ba = sample_data$ba,
  r = c(3, 4),  # Radius within the spatial scale
  dist_weighted = TRUE
)

neigh_multi_r(
  sp = sample_data$latin,
  gx = sample_data$gx,
  gy = sample_data$gy,
  ba = sample_data$ba,
  r_values = c(3, 5, 10)
)
