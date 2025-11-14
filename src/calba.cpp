#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Template function to process a focal tree for basal area calculations
template <typename DecayFunc>
void process_focal_tree_trimmed(int i, const StringVector& sp, const NumericVector& gx,
                                const NumericVector& gy, const NumericVector& ba, double r,
                                DecayFunc decay_func,
                                NumericVector& con_ba, NumericVector& total_ba) {
  String target_sp = sp[i];
  double con_ba_sum = 0.0;
  double total_ba_sum = 0.0;

  double gx_i = gx[i];
  double gy_i = gy[i];
  int n_focal = gx.size();

  for (int j = 0; j < n_focal; j++) {
    if (i == j) continue;

    // Trim based on bounding box
    double dx = gx[j] - gx_i;
    double dy = gy[j] - gy_i;
    if (std::abs(dx) > r || std::abs(dy) > r) continue;

    // Compute actual distance
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist > 0 && dist <= r) {
      double ba_decay = decay_func(ba[j], dist);

      total_ba_sum += ba_decay; // Add to total basal area
      if (sp[j] == target_sp) {
        con_ba_sum += ba_decay;
      }
    }
  }

  con_ba[i] = con_ba_sum;
  total_ba[i] = total_ba_sum;
}

struct RadiusEntry {
  double radius;
  int index;
};

// [[Rcpp::export]]
List calculate_neighborhood_multi_radius(StringVector sp, NumericVector gx, NumericVector gy,
                                         NumericVector ba, NumericVector r_values,
                                         bool dist_weighted = false) {
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

  for (R_xlen_t i = 0; i < n_focal; ++i) {
    String target_sp = sp[i];
    double gx_i = gx[i];
    double gy_i = gy[i];

    for (R_xlen_t j = 0; j < n_focal; ++j) {
      if (i == j) continue;

      double dx = gx[j] - gx_i;
      double dy = gy[j] - gy_i;
      if (std::abs(dx) > max_r || std::abs(dy) > max_r) continue;

      double dist = std::sqrt(dx * dx + dy * dy);
      if (dist <= 0 || dist > max_r) continue;

      double contribution = dist_weighted ? ba[j] / dist : ba[j];
      auto comparator = [](const RadiusEntry& entry, double value) {
        return entry.radius < value;
      };
      auto it = std::lower_bound(radii.begin(), radii.end(), dist, comparator);
      for (auto iter = it; iter != radii.end(); ++iter) {
        update_matrix(i, contribution, sp[j] == target_sp, *iter);
      }
    }
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
                                 NumericVector ba, double r, bool dist_weighted = false) {
  int n_focal = gx.size();
  NumericVector con_ba(n_focal, 0.0);
  NumericVector total_ba(n_focal, 0.0);

  // Decay function for simple decay
  auto decay_func = [dist_weighted](double ba_j, double dist) {
    if (dist_weighted) {
      return ba_j / dist; // Distance-weighted calculation
    } else {
      return ba_j; // Simple basal area
    }
  };

  // Iterate over focal trees
  for (int i = 0; i < n_focal; i++) {
    process_focal_tree_trimmed(i, sp, gx, gy, ba, r, decay_func, con_ba, total_ba);
  }

  return List::create(
    Named("con_ba") = con_ba,
    Named("total_ba") = total_ba
  );
}

// [[Rcpp::export]]
List calculate_basal_area_decay(NumericVector mu_values, StringVector sp, NumericVector gx,
                                NumericVector gy, NumericVector ba, double r, std::string decay_type) {
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

  for (int i = 0; i < n_focal; ++i) {
    double gx_i = gx[i];
    double gy_i = gy[i];
    double ba_i = ba[i];

    for (int j = i + 1; j < n_focal; ++j) {
      double dx = gx[j] - gx_i;
      double dy = gy[j] - gy_i;

      if (std::abs(dx) > r || std::abs(dy) > r) continue;

      double dist_sq = dx * dx + dy * dy;
      if (dist_sq <= 0 || dist_sq > r_sq) continue;

      double dist = std::sqrt(dist_sq);
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
        total_ba_matrix(i, m) += contrib_i;
        if (same_species) {
          con_ba_matrix(i, m) += contrib_i;
        }

        double contrib_j = ba_i * decay_value;
        total_ba_matrix(j, m) += contrib_j;
        if (same_species) {
          con_ba_matrix(j, m) += contrib_j;
        }
      }
    }
  }

  return List::create(
    Named("con_ba_matrix") = con_ba_matrix,
    Named("total_ba_matrix") = total_ba_matrix
  );
}

// [[Rcpp::export]]
NumericVector count_total_cpp(NumericVector gx, NumericVector gy, double r) {
  int n = gx.size();
  NumericVector res(n);

  for (int j = 0; j < n; j++) {
    int trees = 0;
    double gx_j = gx[j];
    double gy_j = gy[j];

    for (int i = 0; i < n; i++) {
      if (i == j) continue; // Skip the focal tree
      double dx = gx[i] - gx_j;
      double dy = gy[i] - gy_j;

      // Trim based on bounding box
      if (std::abs(dx) > r || std::abs(dy) > r) continue;

      // Compute actual distance
      if ((dx * dx + dy * dy) <= r * r) {
        trees++;
      }
    }
    res[j] = trees;
  }

  return res;
}

// [[Rcpp::export]]
NumericVector count_con_cpp(StringVector sp, NumericVector gx, NumericVector gy, double r) {
  int n = sp.size();
  NumericVector res(n);

  for (int j = 0; j < n; j++) {
    int trees = 0;
    double gx_j = gx[j];
    double gy_j = gy[j];
    String target_sp = sp[j];

    for (int i = 0; i < n; i++) {
      if (i == j) continue; // Skip the focal tree
      double dx = gx[i] - gx_j;
      double dy = gy[i] - gy_j;

      // Trim based on bounding box
      if (std::abs(dx) > r || std::abs(dy) > r) continue;

      // Compute actual distance and check species
      if ((dx * dx + dy * dy) <= r * r && sp[i] == target_sp) {
        trees++;
      }
    }
    res[j] = trees;
  }

  return res;
}
