#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
using namespace Rcpp;

struct RadiusEntry {
  double radius;
  int index;
};

struct SpatialGrid {
  double xmin;
  double ymin;
  double xmax;
  double ymax;
  double cell_size;
  int nx;
  int ny;
  std::vector<std::vector<int>> cells;
  std::vector<int> cell_ids;
};

enum class EdgeCorrection {
  None,
  Safe
};

EdgeCorrection parse_edge_correction(const std::string& value) {
  if (value == "none") {
    return EdgeCorrection::None;
  }
  if (value == "safe") {
    return EdgeCorrection::Safe;
  }
  stop("`edge_correction` must be either 'none' or 'safe'");
}

inline bool is_safe_focal(double x, double y, double r,
                          double xmin, double xmax,
                          double ymin, double ymax) {
  return (x - r >= xmin) && (x + r <= xmax) &&
         (y - r >= ymin) && (y + r <= ymax);
}

inline int clamp_value(int value, int lower, int upper) {
  if (value < lower) {
    return lower;
  }
  if (value > upper) {
    return upper;
  }
  return value;
}

inline int compute_cell_id(const SpatialGrid& grid, double x, double y) {
  int cx = static_cast<int>(std::floor((x - grid.xmin) / grid.cell_size));
  int cy = static_cast<int>(std::floor((y - grid.ymin) / grid.cell_size));
  cx = clamp_value(cx, 0, grid.nx - 1);
  cy = clamp_value(cy, 0, grid.ny - 1);
  return cy * grid.nx + cx;
}

SpatialGrid build_spatial_grid(const NumericVector& gx, const NumericVector& gy, double radius) {
  int n = gx.size();
  SpatialGrid grid;
  grid.cell_size = radius;
  if (n == 0) {
    grid.xmin = 0.0;
    grid.ymin = 0.0;
    grid.nx = 1;
    grid.ny = 1;
    grid.cells.assign(1, std::vector<int>());
    grid.cell_ids.clear();
    return grid;
  }
  grid.xmin = gx[0];
  grid.ymin = gy[0];
  double max_x = gx[0];
  double max_y = gy[0];

  for (int i = 1; i < n; ++i) {
    if (gx[i] < grid.xmin) grid.xmin = gx[i];
    if (gx[i] > max_x) max_x = gx[i];
    if (gy[i] < grid.ymin) grid.ymin = gy[i];
    if (gy[i] > max_y) max_y = gy[i];
  }

  grid.xmax = max_x;
  grid.ymax = max_y;

  double range_x = max_x - grid.xmin;
  double range_y = max_y - grid.ymin;
  double max_range = std::max(range_x, range_y);
  double approx_cells = std::sqrt(static_cast<double>(n));
  if (approx_cells < 1.0) {
    approx_cells = 1.0;
  }
  double scale = (max_range > 0.0) ? (max_range / approx_cells) : radius;
  double cell_size = std::max(radius, scale);
  grid.cell_size = cell_size;

  int nx = static_cast<int>(std::floor(range_x / cell_size)) + 1;
  int ny = static_cast<int>(std::floor(range_y / cell_size)) + 1;
  grid.nx = std::max(1, nx);
  grid.ny = std::max(1, ny);

  grid.cells.assign(grid.nx * grid.ny, std::vector<int>());
  grid.cell_ids.assign(n, 0);

  for (int i = 0; i < n; ++i) {
    int id = compute_cell_id(grid, gx[i], gy[i]);
    grid.cell_ids[i] = id;
    grid.cells[id].push_back(i);
  }

  return grid;
}

template <typename Func>
void for_each_neighbor(const SpatialGrid& grid, int focal_idx, const NumericVector& gx,
                       const NumericVector& gy, double max_radius_sq, Func&& func) {
  int base_cell = grid.cell_ids[focal_idx];
  int base_cx = base_cell % grid.nx;
  int base_cy = base_cell / grid.nx;
  double gx_i = gx[focal_idx];
  double gy_i = gy[focal_idx];

  for (int dx = -1; dx <= 1; ++dx) {
    int ncx = base_cx + dx;
    if (ncx < 0 || ncx >= grid.nx) continue;
    for (int dy = -1; dy <= 1; ++dy) {
      int ncy = base_cy + dy;
      if (ncy < 0 || ncy >= grid.ny) continue;
      int cell_index = ncy * grid.nx + ncx;
      const std::vector<int>& bucket = grid.cells[cell_index];
      for (int j : bucket) {
        if (j == focal_idx) continue;
        double dx_val = gx[j] - gx_i;
        double dy_val = gy[j] - gy_i;
        double dist_sq = dx_val * dx_val + dy_val * dy_val;
        if (dist_sq > max_radius_sq) continue;
        double dist = std::sqrt(dist_sq);
        func(j, dist, dist_sq);
      }
    }
  }
}

// [[Rcpp::export]]
List calculate_neighborhood_multi_radius(StringVector sp, NumericVector gx, NumericVector gy,
                                         NumericVector ba, NumericVector r_values,
                                         bool dist_weighted = false,
                                         std::string edge_correction = "none") {
  R_xlen_t n_focal = gx.size();
  R_xlen_t n_r = r_values.size();

  if (n_r == 0) {
    stop("`r_values` must contain at least one radius");
  }

  std::vector<RadiusEntry> radii;
  radii.reserve(static_cast<size_t>(n_r));
  for (R_xlen_t k = 0; k < n_r; ++k) {
    double r = r_values[k];
    if (!std::isfinite(r) || r <= 0) {
      stop("`r_values` must be positive and finite");
    }
    radii.push_back({r, static_cast<int>(k)});
  }

  std::sort(radii.begin(), radii.end(), [](const RadiusEntry& a, const RadiusEntry& b) {
    return a.radius < b.radius;
  });

  double max_r = radii.back().radius;

  NumericMatrix con_ba(n_focal, n_r);
  NumericMatrix total_ba(n_focal, n_r);
  NumericMatrix con_count(n_focal, n_r);
  NumericMatrix total_count(n_focal, n_r);

  auto update_matrix = [&](R_xlen_t idx, double contribution, bool same_species, const RadiusEntry& entry) {
    int col = entry.index;
    total_ba(idx, col) += contribution;
    total_count(idx, col) += 1;
    if (same_species) {
      con_ba(idx, col) += contribution;
      con_count(idx, col) += 1;
    }
  };

  SpatialGrid grid = build_spatial_grid(gx, gy, max_r);
  double max_r_sq = max_r * max_r;

  EdgeCorrection correction = parse_edge_correction(edge_correction);
  bool require_safe = (correction == EdgeCorrection::Safe);
  double na_value = NA_REAL;

  auto comparator = [](const RadiusEntry& entry, double value) {
    return entry.radius < value;
  };

  for (R_xlen_t i = 0; i < n_focal; ++i) {
    String target_sp = sp[i];
    bool safe = is_safe_focal(gx[i], gy[i], max_r, grid.xmin, grid.xmax, grid.ymin, grid.ymax);
    if (require_safe && !safe) {
      for (int col = 0; col < n_r; ++col) {
        con_ba(i, col) = na_value;
        total_ba(i, col) = na_value;
        con_count(i, col) = na_value;
        total_count(i, col) = na_value;
      }
      continue;
    }
    for_each_neighbor(grid, i, gx, gy, max_r_sq, [&](int j, double dist, double /*dist_sq*/) {
      if (dist <= 0) {
        return;
      }
      double contribution = dist_weighted ? ba[j] / dist : ba[j];
      auto it = std::lower_bound(radii.begin(), radii.end(), dist, comparator);
      for (auto iter = it; iter != radii.end(); ++iter) {
        update_matrix(i, contribution, sp[j] == target_sp, *iter);
      }
    });
  }

  return List::create(
    Named("r_values") = r_values,
    Named("con_ba") = con_ba,
    Named("total_ba") = total_ba,
    Named("con_count") = con_count,
    Named("total_count") = total_count
  );
}

// [[Rcpp::export]]
List calculate_basal_area_simple(StringVector sp, NumericVector gx, NumericVector gy,
                                 NumericVector ba, double r, bool dist_weighted = false,
                                 std::string edge_correction = "none") {
  int n_focal = gx.size();
  NumericVector con_ba(n_focal, 0.0);
  NumericVector total_ba(n_focal, 0.0);

  SpatialGrid grid = build_spatial_grid(gx, gy, r);
  double max_radius_sq = r * r;

  EdgeCorrection correction = parse_edge_correction(edge_correction);
  bool require_safe = (correction == EdgeCorrection::Safe);
  double na_value = NA_REAL;

  for (int i = 0; i < n_focal; ++i) {
    bool safe = is_safe_focal(gx[i], gy[i], r, grid.xmin, grid.xmax, grid.ymin, grid.ymax);
    if (require_safe && !safe) {
      con_ba[i] = na_value;
      total_ba[i] = na_value;
      continue;
    }
    String target_sp = sp[i];
    for_each_neighbor(grid, i, gx, gy, max_radius_sq, [&](int j, double dist, double /*dist_sq*/) {
      if (dist <= 0) {
        return;
      }
      double contribution = dist_weighted ? ba[j] / dist : ba[j];
      total_ba[i] += contribution;
      if (sp[j] == target_sp) {
        con_ba[i] += contribution;
      }
    });
  }

  return List::create(
    Named("con_ba") = con_ba,
    Named("total_ba") = total_ba
  );
}

// [[Rcpp::export]]
List calculate_basal_area_decay(NumericVector mu_values, StringVector sp, NumericVector gx,
                                NumericVector gy, NumericVector ba, double r, std::string decay_type,
                                std::string edge_correction = "none") {
  int n_focal = gx.size();
  int n_mu = mu_values.size();
  NumericMatrix con_ba_matrix(n_focal, n_mu);
  NumericMatrix total_ba_matrix(n_focal, n_mu);

  if (n_mu == 0) {
    stop("`mu_values` must contain at least one value");
  }

  bool exponential_decay = (decay_type == "exponential");
  bool exponential_normal_decay = (decay_type == "exponential-normal");
  if (!exponential_decay && !exponential_normal_decay) {
    stop("Unknown decay type");
  }

  std::vector<double> inv_mu(n_mu);
  std::vector<double> inv_mu_sq(n_mu);
  for (int m = 0; m < n_mu; ++m) {
    double mu = mu_values[m];
    if (mu <= 0 || !std::isfinite(mu)) {
      stop("`mu_values` must be positive and finite");
    }
    inv_mu[m] = exponential_decay ? 1.0 / mu : 0.0;
    inv_mu_sq[m] = exponential_normal_decay ? 1.0 / (mu * mu) : 0.0;
  }

  double r_sq = r * r;
  SpatialGrid grid = build_spatial_grid(gx, gy, r);
  EdgeCorrection correction = parse_edge_correction(edge_correction);
  bool require_safe = (correction == EdgeCorrection::Safe);
  double na_value = NA_REAL;

  for (int i = 0; i < n_focal; ++i) {
    double ba_i = ba[i];
    bool safe_i = is_safe_focal(gx[i], gy[i], r, grid.xmin, grid.xmax, grid.ymin, grid.ymax);
    if (require_safe && !safe_i) {
      for (int m = 0; m < n_mu; ++m) {
        con_ba_matrix(i, m) = na_value;
        total_ba_matrix(i, m) = na_value;
      }
    }
    for_each_neighbor(grid, i, gx, gy, r_sq, [&](int j, double dist, double dist_sq) {
      if (j <= i || dist_sq <= 0) {
        return;
      }
      double ba_j = ba[j];
      bool same_species = sp[j] == sp[i];
      for (int m = 0; m < n_mu; ++m) {
        double decay_value;
        if (exponential_decay) {
          decay_value = std::exp(-dist * inv_mu[m]);
        } else {
          decay_value = std::exp(-dist_sq * inv_mu_sq[m]);
        }

        double contrib_i = ba_j * decay_value;
        if (!(require_safe && !safe_i)) {
          total_ba_matrix(i, m) += contrib_i;
          if (same_species) {
            con_ba_matrix(i, m) += contrib_i;
          }
        }

        double contrib_j = ba_i * decay_value;
        total_ba_matrix(j, m) += contrib_j;
        if (same_species) {
          con_ba_matrix(j, m) += contrib_j;
        }
      }
    });
  }

  return List::create(
    Named("con_ba_matrix") = con_ba_matrix,
    Named("total_ba_matrix") = total_ba_matrix
  );
}

// [[Rcpp::export]]
NumericVector count_total_cpp(NumericVector gx, NumericVector gy, double r,
                             std::string edge_correction = "none") {
  int n = gx.size();
  NumericVector res(n);

  SpatialGrid grid = build_spatial_grid(gx, gy, r);
  double max_radius_sq = r * r;
  EdgeCorrection correction = parse_edge_correction(edge_correction);
  bool require_safe = (correction == EdgeCorrection::Safe);
  double na_value = NA_REAL;

  for (int j = 0; j < n; ++j) {
    if (require_safe) {
      bool safe = is_safe_focal(gx[j], gy[j], r, grid.xmin, grid.xmax, grid.ymin, grid.ymax);
      if (!safe) {
        res[j] = na_value;
        continue;
      }
    }
    int trees = 0;
    for_each_neighbor(grid, j, gx, gy, max_radius_sq, [&](int /*i*/, double /*dist*/, double dist_sq) {
      if (dist_sq <= 0.0) {
        return;
      }
      trees++;
    });
    res[j] = trees;
  }

  return res;
}

// [[Rcpp::export]]
NumericVector count_con_cpp(StringVector sp, NumericVector gx, NumericVector gy, double r,
                           std::string edge_correction = "none") {
  int n = sp.size();
  NumericVector res(n);

  SpatialGrid grid = build_spatial_grid(gx, gy, r);
  double max_radius_sq = r * r;
  EdgeCorrection correction = parse_edge_correction(edge_correction);
  bool require_safe = (correction == EdgeCorrection::Safe);
  double na_value = NA_REAL;

  for (int j = 0; j < n; ++j) {
    int trees = 0;
    String target_sp = sp[j];
    if (require_safe) {
      bool safe = is_safe_focal(gx[j], gy[j], r, grid.xmin, grid.xmax, grid.ymin, grid.ymax);
      if (!safe) {
        res[j] = na_value;
        continue;
      }
    }
    for_each_neighbor(grid, j, gx, gy, max_radius_sq, [&](int i, double /*dist*/, double dist_sq) {
      if (dist_sq <= 0.0) {
        return;
      }
      if (sp[i] == target_sp) {
        trees++;
      }
    });
    res[j] = trees;
  }

  return res;
}
