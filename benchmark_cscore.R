library(microbenchmark)
library(stats)

# Mock function from package
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
    flanking_values <- c(rep(NA, flanking), flanking_values)
    return(flanking_values)
}

# Proposed zoo implementation
compute_vec_cscore_zoo <- function(count, flanking = 6, method) {
    # We need to exclude the center.
    # Window size = 2 * flanking + 1
    # But we want excludes center.
    # For mean: sum(window) - center / (2*flanking)
    # For median: rolling median of '2*flanking' elements excluding center.
    # zoo::rollmedian doesn't directly support excluding center efficiently without tricks.

    n <- length(count)
    res <- rep(NA, n)

    if (method == "mean") {
        # rolling sum of 2*flanking+1
        # subtract center
        # divide by 2*flanking
        # zoo::rollsum
        k <- 2 * flanking + 1
        # fill=NA to handle edges
        sums <- zoo::rollsum(count, k = k, fill = NA, align = "center")
        # center values at indices where sums is not NA
        # indices valid: (flanking+1) to (n-flanking)
        center_idx <- (flanking + 1):(n - flanking)

        # sums[center_idx] includes count[center_idx].
        # We remove it.
        means <- (sums[center_idx] - count[center_idx]) / (2 * flanking)
        res[center_idx] <- means
    } else if (method == "median") {
        # Excluding center for median is trickier with standard rolling functions.
        # But median is robust. Maybe including center doesn't change much?
        # User requirement: "The C-score represents a drop ... compared to the environmental coverage"
        # Excluding center is explicit in current code.

        # If we MUST exclude center, zoo is harder for median.
        # We can try to see performance of just 'apply' over matrix of windows (embed)
        # embed is also memory heavy.

        # If zoo cannot do it easily, Rcpp is the way.
        pass
    }
    return(res)
}

# Rcpp candidate prototype (if we were to use it)
# We won't benchmark Rcpp implementation here as we don't want to compile yet.

# Data
set.seed(42)
N <- 100000
count <- runif(N, 0, 100)

cat("Loading package functions...\n")
devtools::load_all(".")
# Now .compute_vec_cscore uses Rcpp

cat("Benchmarking Mean (N=100,000)...\n")
try({
    print(microbenchmark(
        orig = compute_vec_cscore_orig(count, 6, "mean"),
        rcpp = .compute_vec_cscore(count, 6, "mean"),
        times = 20
    ))
})

cat("\nBenchmarking Median (N=100,000)...\n")
try({
    print(microbenchmark(
        orig = compute_vec_cscore_orig(count, 6, "median"),
        rcpp = .compute_vec_cscore(count, 6, "median"),
        times = 20
    ))
})
