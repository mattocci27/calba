library(tidyverse)
library(hexSticker)
library(showtext)
font_add_google("Nunito", "nunito")
showtext_auto()

filter_hex <- function(
  df, x_col, y_col,
  s_x = 1.0, s_y = 0.9, s_width = 1.55, s_height = 1.55,
  hex_cx = 1, hex_cy = 1, hex_r = 1.0, margin = 0.01
) {
  point_in_poly <- function(px, py, vx, vy) {
    n <- length(vx)
    inside <- logical(length(px))
    j <- n
    for (i in seq_len(n)) {
      xi <- vx[i]; yi <- vy[i]
      xj <- vx[j]; yj <- vy[j]
      intersect <- ((yi > py) != (yj > py)) &
        (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
      inside <- xor(inside, intersect)
      j <- i
    }
    inside
  }

  # pointy-top hexagon in sticker coordinates; small margin to sit safely inside the outline
  theta <- seq(0, 2 * pi, length.out = 7)[1:6] + pi / 2
  hex_r_eff <- max(hex_r - margin, 0)
  vx_hex <- hex_cx + hex_r_eff * cos(theta)
  vy_hex <- hex_cy + hex_r_eff * sin(theta)

  # transform to plot coordinates (data assumed on [0,1] x [0,1])
  vx_plot <- (vx_hex - s_x) / s_width + 0.5
  vy_plot <- (vy_hex - s_y) / s_height + 0.5

  inside <- point_in_poly(df[[x_col]], df[[y_col]], vx_plot, vy_plot)
  df[inside, , drop = FALSE]
}

set.seed(123)

# ----- focal tree and radii -----

focal   <- data.frame(x = 0.30, y = 0.30)
r_small <- 0.18   # focal circle radius

# share subplot placement between filter and sticker
s_x      <- 1.0
s_y      <- 0.9
s_width  <- 1.55
s_height <- 1.55
h_size   <- 1.0   # hex radius used by both filter and sticker

# ----- generate many candidate points, then cluster to 12 centres -----

n_candidates <- 200

candidates <- data.frame(
  x = runif(n_candidates, 0, 1),
  y = runif(n_candidates, 0, 1)
)

km <- kmeans(candidates[, c("x", "y")], centers = 12, iter.max = 50)
trees <- as.data.frame(km$centers)

# ----- avoid clutter inside the small focal circle -----

dx <- trees$x - focal$x
dy <- trees$y - focal$y
dist <- sqrt(dx^2 + dy^2)

inside <- dist < r_small * 1.1

# move points that are too close to the focal tree out to just outside the circle
trees[inside, c("x", "y")] <- data.frame(
  x = focal$x + dx[inside] / dist[inside] * r_small * 1.1,
  y = focal$y + dy[inside] / dist[inside] * r_small * 1.1
)


# ----- assign species labels in a balanced way (6 and 6) -----

trees$sp <- rep(c("A", "B"), length.out = nrow(trees))

# ----- clip to the actual hex outline so nothing leaks out -----

trees <- filter_hex(
  df      = trees,
  x_col   = "x",
  y_col   = "y",
  s_x     = s_x,
  s_y     = s_y,
  s_width = s_width,
  s_height = s_height,
  hex_r   = h_size
)

print(trees)

trees <- trees |>
  mutate(y = ifelse(y > 0.8, y + 0.05, y)) |>
  mutate(x = ifelse(row_number() == 3, x + 0.015, x))

# ----- colour palette -----

col_bg_hex <- "#f9fafc"   # hex background
col_border <- "#203040"   # hex border and package name
col_sp_A   <- "#2c7f5e"   # species A (green)
col_sp_B   <- "#3b6fb6"   # species B (blue)
col_focal  <- "#f39c12"   # focal tree (orange)
col_circle <- "#f5b041"   # radius circle (lighter orange)

# ----- ggplot for the inner artwork -----

theta <- pi / 3   # choose angle of the radius (here 60 degrees)

p <- ggplot(trees, aes(x, y)) +
  # radius line from center to circle
  annotate(
    "segment",
    x = focal$x,
    y = focal$y,
    xend = focal$x + r_small * cos(theta),
    yend = focal$y + r_small * sin(theta),
    color = col_circle,
    linewidth = 0.8,
    alpha = 0.9
  ) +
  geom_point(
    aes(color = sp),
    size = 2.4,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(A = col_sp_A, B = col_sp_B)) +
  geom_point(
    data  = focal,
    size  = 2.6,
    color = col_focal,
    alpha = 0.95
  ) +
  annotate(
    "path",
    x = focal$x + r_small * cos(seq(0, 2 * pi, length.out = 201)),
    y = focal$y + r_small * sin(seq(0, 2 * pi, length.out = 201)),
    color     = col_circle,
    linewidth = 0.7,
    alpha     = 0.96
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = 0) +
  theme_void()


# ----- hex sticker output -----
sticker(
  subplot   = p,
  package   = "calba",
  p_size    = 54,
  p_y       = 1.4,
  p_color   = col_border,
  p_family  = "nunito",
  p_fontface = "bold",
  s_x       = s_x,
  s_y       = s_y,
  s_width   = s_width,
  s_height  = s_height,
  h_size    = h_size,
  h_fill    = col_bg_hex,
  h_color   = col_border,
  white_around_sticker = FALSE,
  dpi       = 600,
  filename  = "man/figures/icon.png"
)

