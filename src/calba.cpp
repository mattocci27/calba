#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Function to precompute the distance matrix
NumericMatrix compute_distance_matrix(const NumericVector& gx, const NumericVector& gy, int n_focal) {
  NumericMatrix dist_matrix(n_focal, n_focal);
  for (int i = 0; i < n_focal; i++) {
    for (int j = 0; j < n_focal; j++) {
      if (i == j) {
        dist_matrix(i, j) = 0; // Distance to self is 0
      } else {
        dist_matrix(i, j) = sqrt(pow(gx[j] - gx[i], 2) + pow(gy[j] - gy[i], 2));
      }
    }
  }
  return dist_matrix;
}

// General function to process a single focal tree
template <typename DecayFunc>
void process_focal_tree(int i, const NumericMatrix& dist_matrix, const CharacterVector& latin,
                        const NumericVector& ba, double r, DecayFunc decay_func,
                        NumericVector& con_ba, NumericVector& het_ba) {
  String target_sp = latin[i];
  double con_ba_sum = 0.0;
  double het_ba_sum = 0.0;

  int n_focal = dist_matrix.ncol();
  for (int j = 0; j < n_focal; j++) {
    if (i == j) continue; // Skip the focal tree itself

    double dist = dist_matrix(i, j);
    if (dist > 0 && dist <= r) {
      double ba_decay = decay_func(ba[j], dist);

      if (latin[j] == target_sp) {
        con_ba_sum += ba_decay;
      } else {
        het_ba_sum += ba_decay;
      }
    }
  }

  con_ba[i] = con_ba_sum;
  het_ba[i] = het_ba_sum;
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
  NumericMatrix het_ba_matrix(n_focal, n_mu);

  // Precompute distances
  NumericMatrix dist_matrix = compute_distance_matrix(gx, gy, n_focal);

  // Iterate over mu values
  for (int m = 0; m < n_mu; m++) {
    double mu = mu_values[m];

    // Decay function for exponential decay
    auto decay_func = [mu](double ba_j, double dist) {
      return ba_j * exp(-dist / mu);
    };

    // Temporary vectors to store results for the current mu
    NumericVector con_ba(n_focal, 0.0);
    NumericVector het_ba(n_focal, 0.0);

    // Iterate over focal trees
    for (int i = 0; i < n_focal; i++) {
      process_focal_tree(i, dist_matrix, latin, ba, r, decay_func, con_ba, het_ba);
    }

    // Store results for this mu
    con_ba_matrix(_, m) = con_ba;
    het_ba_matrix(_, m) = het_ba;
  }

  return List::create(
    Named("con_ba_matrix") = con_ba_matrix,
    Named("het_ba_matrix") = het_ba_matrix
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
  NumericVector het_ba(n_focal, 0.0);

  // Precompute distances
  NumericMatrix dist_matrix = compute_distance_matrix(gx, gy, n_focal);

  // Decay function for simple decay (ba / dist)
  auto decay_func = [](double ba_j, double dist) {
    return ba_j / dist;
  };

  // Iterate over focal trees
  for (int i = 0; i < n_focal; i++) {
    process_focal_tree(i, dist_matrix, latin, ba, r, decay_func, con_ba, het_ba);
  }

  return List::create(
    Named("con_ba") = con_ba,
    Named("het_ba") = het_ba
  );
}

