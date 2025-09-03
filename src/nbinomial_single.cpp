#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

// Optimized helper function to calculate p-value and effect size
inline void calculateFlankStats(double center_reads, double center_pred_bias,
                                double flank_reads, double flank_pred_bias,
                                double& pval, double& effect_size) {
  if (flank_reads == 0) {
    pval = 1.0;
    effect_size = 0.0;
  } else {
    double bias_ratio = flank_pred_bias / (flank_pred_bias + center_pred_bias);
    pval = R::pnbinom(center_reads, flank_reads, bias_ratio, 1, 0);
    effect_size = -std::log((center_reads + 1e-10) * (flank_pred_bias + center_pred_bias + 1e-10) /
      ((flank_reads + center_reads + 1e-10) * center_pred_bias));
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix footprintingAnalysis(NumericVector center_sizes, int center_position, int start, int end,
                                   NumericVector flank_windows, bool wide_flank, DataFrame bardata) {
  int n = center_sizes.size();
  NumericMatrix results(n, 2); // Stores pval and effect_size for each center_size

  // Extract columns from bardata
  IntegerVector positions = bardata["position"];
  NumericVector Tn5Insertion = bardata["Tn5Insertion"];
  NumericVector pred_bias = bardata["pred_bias"];

  // Build a position-to-index lookup table
  std::unordered_map<int, int> position_map;
  for (int i = 0; i < positions.size(); ++i) {
    position_map[positions[i]] = i;
  }

  // Load ACAT function from R
  Rcpp::Function ACAT("ACAT", Rcpp::Environment::namespace_env("ACAT"));

  for (int i = 0; i < n; ++i) {
    double center_size = center_sizes[i];
    double pval = 1.0, effect_size = 0.0;

    // Skip calculations if too close to start/end
    if (std::abs(center_position - start) <= center_size || std::abs(end - center_position) <= center_size) {
      results(i, 0) = pval;
      results(i, 1) = effect_size;
      continue;
    }

    // Compute flank sizes (avoiding pmin and pmax for performance)
    NumericVector flank_sizes = center_size + flank_windows;
    // NumericVector flank_sizes = 1;
    // NumericVector flank_sizes = 20;
    // flank_sizes[flank_sizes > 100] = 100;
    flank_sizes[flank_sizes < 1] = 1;
    flank_sizes = unique(flank_sizes);
    if (wide_flank) {
      flank_sizes = 2 * flank_sizes;
    }

    // Compute Center Region
    double center_reads = 0.0, center_pred_bias = 0.0;
    for (int j = center_position - center_size; j <= center_position + center_size; ++j) {
      auto it = position_map.find(j);
      if (it != position_map.end()) {
        int idx = it->second;
        center_reads += Tn5Insertion[idx];
        center_pred_bias += pred_bias[idx];
      }
    }

    NumericVector pvals(flank_sizes.size()), effect_sizes(flank_sizes.size());

    // Process each flank_size
    for (int j = 0; j < flank_sizes.size(); ++j) {
      double flank_size = flank_sizes[j];

      // Calculate Flank Boundaries (avoid redundant min/max computations)
      int left_start = std::max(start, static_cast<int>(center_position - center_size - flank_size));
      int left_end = std::max(start, static_cast<int>(center_position - center_size - 1));
      int right_start = std::min(end, static_cast<int>(center_position + center_size + 1));
      int right_end = std::min(end, static_cast<int>(center_position + center_size + flank_size));

      // Compute Left & Right Flanks
      double left_reads = 0.0, right_reads = 0.0;
      double left_pred_bias = 0.0, right_pred_bias = 0.0;

      for (int k = left_start; k <= left_end; ++k) {
        auto it = position_map.find(k);
        if (it != position_map.end()) {
          int idx = it->second;
          left_reads += Tn5Insertion[idx];
          left_pred_bias += pred_bias[idx];
        }
      }

      for (int k = right_start; k <= right_end; ++k) {
        auto it = position_map.find(k);
        if (it != position_map.end()) {
          int idx = it->second;
          right_reads += Tn5Insertion[idx];
          right_pred_bias += pred_bias[idx];
        }
      }

      // Calculate p-values and effect sizes
      double pval_left, pval_right, effect_size_left, effect_size_right;
      calculateFlankStats(center_reads, center_pred_bias, left_reads, left_pred_bias, pval_left, effect_size_left);
      calculateFlankStats(center_reads, center_pred_bias, right_reads, right_pred_bias, pval_right, effect_size_right);

      pvals[j] = std::max(std::max(pval_left, pval_right), 1e-50);
      effect_sizes[j] = (pval_left > pval_right) ? effect_size_left : effect_size_right;

    }

    // Adjust pvals and compute weights
    pvals[pvals == 1.0] = 0.95;  // Avoid perfect p-values
    NumericVector weights = exp(-abs(flank_sizes - center_size));
    weights = weights / sum(weights); // Normalize

    // Use ACAT to combine p-values
    pval = Rcpp::as<double>(ACAT(pvals, Named("weights") = weights));
    effect_size = sum(effect_sizes * weights);
    results(i, 0) = pval;
    results(i, 1) = effect_size;
  }

  return results;
}

