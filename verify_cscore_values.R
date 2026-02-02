# Script to verify exactness of C-score computation results
# Compares original R implementation vs New Rcpp implementation

library(testthat)
library(stats)

# 1. Define Original R Implementation
compute_vec_cscore_orig <- function(count, flanking = 6, method) {
    if (method == "median") {
        flanking_values <- vapply(seq(1 + flanking, length(count)), function(x) {
            median(count[c((x - flanking):(x - 1), (x + 1):(x + flanking))])
        }, numeric(1))
    } else if (method == "mean") {
        flanking_values <- vapply(seq(1 + flanking, length(count)), function(x) {
            mean(count[c((x - flanking):(x - 1), (x + 1):(x + flanking))])
        }, numeric(1))
    }
    # Prepend NA
    flanking_values <- c(rep(NA, flanking), flanking_values)

    if (method %in% c("median", "mean")) {
        scorec_raw <- 1 - count / flanking_values
    }
    return(pmax(scorec_raw, 0))
}

# 2. Load New Implementation (Rcpp)
devtools::load_all(".", quiet = TRUE)
# The package function .compute_vec_cscore now uses Rcpp

# 3. Generate Test Data
set.seed(123)
N <- 1000 # Test vector size
count_data <- runif(N, 0, 1000)
# Add some zeros and potential edge cases
count_data[10:20] <- 0
count_data[50] <- 1000

# 4. Compare Results for MEAN
cat("Comparing MEAN method...\n")
res_orig_mean <- compute_vec_cscore_orig(count_data, flanking = 6, method = "mean")
res_new_mean <- .compute_vec_cscore(count_data, flanking = 6, method = "mean")

# Check equality (ignoring NAs at start)
# We expect NAs at the first 6 positions, and potentially at end if window goes out of bounds?
# The original code: seq(1+flanking, length(count)).
# If length(count) < 1+flanking loop doesn't run.
# Original code logic:
# i goes from 7 to 1000.
# window [i-6, i-1] and [i+1, i+6].
# At i=1000, i+6 = 1006 -> out of bounds?
# R indexing: count[1001:1006] return NAs. mean(..., na.rm=FALSE) returns NA.
# So original implementation returns NA for last 6 positions too?
# Let's check exactness.

equal_mean <- all.equal(res_orig_mean, res_new_mean, tolerance = 1e-10)
if (isTRUE(equal_mean)) {
    cat("[PASS] MEAN results are identical.\n")
} else {
    cat("[FAIL] MEAN results differ!\n")
    print(equal_mean)
    # Show diff
    diff_idx <- which(abs(res_orig_mean - res_new_mean) > 1e-10)
    head(diff_idx)
}

# 5. Compare Results for MEDIAN
cat("\nComparing MEDIAN method...\n")
res_orig_med <- compute_vec_cscore_orig(count_data, flanking = 6, method = "median")
res_new_med <- .compute_vec_cscore(count_data, flanking = 6, method = "median")

equal_med <- all.equal(res_orig_med, res_new_med, tolerance = 1e-10)
if (isTRUE(equal_med)) {
    cat("[PASS] MEDIAN results are identical.\n")
} else {
    cat("[FAIL] MEDIAN results differ!\n")
    print(equal_med)
}

# 6. Edge Case: Small Vector
cat("\nComparing Small Vector (N=10)...\n")
small_count <- runif(10)
# Flanking 6 -> 2*6+1 = 13 window size needed?
# Original loop: 1+6 = 7 to 10.
# i=7: window [1-6] -> [1,6]. [8-13] -> [8,10] + NAs.
# median with NAs returns NA (default).
# Rcpp implementation handles bounds explicitly?
res_small_orig <- compute_vec_cscore_orig(small_count, flanking = 6, method = "median")
res_small_new <- .compute_vec_cscore(small_count, flanking = 6, method = "median")

# Both should likely be all NAs or mostly NAs
print(rbind(original = res_small_orig, new = res_small_new))

equal_small <- all.equal(res_small_orig, res_small_new)
if (isTRUE(equal_small)) {
    cat("[PASS] Small vector results match.\n")
} else {
    cat("[FAIL] Small vector results differ.\n")
}

cat("\nVerification Complete.\n")
