#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Template function to process a focal tree for basal area calculations
template <typename DecayFunc>
void process_focal_tree_trimmed(int i, const StringVector& sp, const NumericVector& gx,
                                const NumericVector& gy, const NumericVector& ba, double r,
                                DecayFunc decay_func, bool dist_weighted,
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
      double ba_decay;

      if (dist_weighted) {
        ba_decay = ba[j] / dist; // Use distance-weighted calculation
      } else {
        ba_decay = decay_func(ba[j], dist); // Use decay function
      }

      total_ba_sum += ba_decay; // Add to total basal area
      if (sp[j] == target_sp) {
        con_ba_sum += ba_decay;
      }
    }
  }

  con_ba[i] = con_ba_sum;
  total_ba[i] = total_ba_sum;
}

// [[Rcpp::export]]
List calculate_basal_area_simple(StringVector sp, NumericVector gx, NumericVector gy,
                                 NumericVector ba, double r, bool dist_weighted = false) {
  int n_focal = gx.size();
  NumericVector con_ba(n_focal, 0.0);
  NumericVector total_ba(n_focal, 0.0);

  // Decay function for simple decay (ba / dist)
  auto decay_func = [](double ba_j, double dist) {
    return ba_j; // This decay_func is effectively bypassed for distance weighting
  };

  // Iterate over focal trees
  for (int i = 0; i < n_focal; i++) {
    process_focal_tree_trimmed(i, sp, gx, gy, ba, r, decay_func, dist_weighted, con_ba, total_ba);
  }

  return List::create(
    Named("con_ba") = con_ba,
    Named("total_ba") = total_ba
  );
}

// [[Rcpp::export]]
List calculate_basal_area_decay(NumericVector mu_values, StringVector sp, NumericVector gx,
                                NumericVector gy, NumericVector ba, double r) {
  int n_focal = gx.size();
  int n_mu = mu_values.size();
  NumericMatrix con_ba_matrix(n_focal, n_mu);
  NumericMatrix total_ba_matrix(n_focal, n_mu);

  for (int m = 0; m < n_mu; m++) {
    double mu = mu_values[m];

    // Decay function for exponential decay
    auto decay_func = [mu](double ba_j, double dist) {
      return ba_j * exp(-dist / mu);
    };

    // Temporary vectors to store results for the current mu
    NumericVector con_ba(n_focal, 0.0);
    NumericVector total_ba(n_focal, 0.0);

    // Iterate over focal trees
    for (int i = 0; i < n_focal; i++) {
      process_focal_tree_trimmed(i, sp, gx, gy, ba, r, decay_func, false, con_ba, total_ba);
    }

    // Store results for this mu
    con_ba_matrix(_, m) = con_ba;
    total_ba_matrix(_, m) = total_ba;
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


