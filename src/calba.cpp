#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Function to compute distances for trees within radius r
template <typename DecayFunc>
void process_focal_tree_trimmed(int i, const CharacterVector& latin, const NumericVector& gx,
                                const NumericVector& gy, const NumericVector& ba, double r,
                                DecayFunc decay_func, NumericVector& con_ba, NumericVector& total_ba) {
  String target_sp = latin[i];
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
      if (latin[j] == target_sp) {
        con_ba_sum += ba_decay;
      }
    }
  }

  con_ba[i] = con_ba_sum;
  total_ba[i] = total_ba_sum;
}

// [[Rcpp::export]]
List calculate_basal_area_decay(NumericVector mu_values, DataFrame data, int n_focal, double r) {
  // Extract columns from data
  CharacterVector latin = data["latin"];
  NumericVector gx = data["gx"];
  NumericVector gy = data["gy"];
  NumericVector ba = data["ba"];

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
      process_focal_tree_trimmed(i, latin, gx, gy, ba, r, decay_func, con_ba, total_ba);
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
List calculate_basal_area_simple(DataFrame data, int n_focal, double r) {
  // Extract columns from data
  CharacterVector latin = data["latin"];
  NumericVector gx = data["gx"];
  NumericVector gy = data["gy"];
  NumericVector ba = data["ba"];

  NumericVector con_ba(n_focal, 0.0);
  NumericVector total_ba(n_focal, 0.0);

  // Decay function for simple decay (ba / dist)
  auto decay_func = [](double ba_j, double dist) {
    return ba_j / dist;
  };

  // Iterate over focal trees
  for (int i = 0; i < n_focal; i++) {
    process_focal_tree_trimmed(i, latin, gx, gy, ba, r, decay_func, con_ba, total_ba);
  }

  return List::create(
    Named("con_ba") = con_ba,
    Named("total_ba") = total_ba
  );
}
