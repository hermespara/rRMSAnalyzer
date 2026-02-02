#include <Rcpp.h>
using namespace Rcpp;

//' Compute C-score window values (Optimized)
//' 
//' @param count Vector of counts
//' @param flanking Size of the flanking window
//' @param method "mean" or "median"
//' @return Vector of flanking values (same length as count)
//' @export
// [[Rcpp::export]]
NumericVector compute_window_score_cpp(NumericVector count, int flanking, String method) {
  int n = count.size();
  NumericVector result(n);
  
  // Initialize with NA
  std::fill(result.begin(), result.end(), NA_REAL);
  
  // Indices loop from (flanking) to (n - flanking - 1)
  // In R 1-based: (1+flanking) to (n-flanking). So indices [flanking, n-flanking-1]
  
  if (n < 2 * flanking + 1) {
    return result; 
  }

  if (method == "mean") {
    for (int i = flanking; i < n - flanking; i++) {
      double sum = 0;
      // sum left window: [i - flanking, i - 1]
      for (int k = i - flanking; k < i; k++) {
        sum += count[k];
      }
      // sum right window: [i + 1, i + flanking]
      for (int k = i + 1; k <= i + flanking; k++) {
        sum += count[k];
      }
      result[i] = sum / (2.0 * flanking);
    }
  } else if (method == "median") {
    int window_size = 2 * flanking;
    std::vector<double> window_vals(window_size);
    
    for (int i = flanking; i < n - flanking; i++) {
      int idx = 0;
      // left
      for (int k = i - flanking; k < i; k++) {
        window_vals[idx++] = count[k];
      }
      // right
      for (int k = i + 1; k <= i + flanking; k++) {
        window_vals[idx++] = count[k];
      }
      
      // Compute median
      // Using nth_element
      // mid point for size 2*flanking (even number)
      // std::median usually: mean of (n/2 - 1) and (n/2) for even size
      // R's median: "For an even number of data points, the median is the mean of the two middle values."
      
      std::nth_element(window_vals.begin(), window_vals.begin() + window_size / 2, window_vals.end());
      double m2 = window_vals[window_size / 2];
      
      std::nth_element(window_vals.begin(), window_vals.begin() + window_size / 2 - 1, window_vals.end());
      double m1 = window_vals[window_size / 2 - 1];
      
      result[i] = (m1 + m2) / 2.0; 
    }
  }
  
  return result;
}
