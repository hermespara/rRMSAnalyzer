# Script to verify C-score computation on ribo_toy dataset

library(SummarizedExperiment)
devtools::load_all(".", quiet = TRUE)

# 1. Define Original R Implementation (Internal Logic)
# We need to replicate the exact logic of .compute_sample_cscore_se but using the slow vapply
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

    if (method %in% c("median", "mean")) {
        scorec_raw <- 1 - count / flanking_values
    }
    return(pmax(scorec_raw, 0))
}

compute_sample_cscore_orig <- function(sample_counts, rna_factor, flanking, method) {
    counts_by_rna <- split(sample_counts, rna_factor)
    scores_by_rna <- lapply(counts_by_rna, compute_vec_cscore_orig, flanking, method)
    scores <- unsplit(scores_by_rna, rna_factor)
    return(scores)
}

# 2. Load Data
data("ribo_toy")
# ribo_toy already has cscores computed? No, usually it has counts.
# Let's recompute from scratch using counts.

counts <- assay(ribo_toy, "counts")
rnas <- rowData(ribo_toy)$rna

cat("Dataset dimensions:", nrow(counts), "sites x", ncol(counts), "samples\n")

# 3. Compute using Original Method (Manual Loop)
cat("Computing C-scores with Original Method (Slow)...\n")
t1 <- Sys.time()
cscores_orig <- apply(counts, 2, compute_sample_cscore_orig, rna_factor = rnas, flanking = 6, method = "median")
t2 <- Sys.time()
cat("Original computation took:", format(t2 - t1), "\n")

# 4. Compute using New Package Method (Rcpp)
cat("Computing C-scores with New Method (Rcpp)...\n")
# We use the exported function compute_cscore which uses the optimized internal function
t3 <- Sys.time()
ribo_new <- compute_cscore(ribo_toy, flanking = 6, method = "median", ncores = 1)
cscores_new <- assay(ribo_new, "cscore")
t4 <- Sys.time()
cat("New computation took:", format(t4 - t3), "\n")

# 5. Compare
cat("\nComparing results...\n")
equal_res <- all.equal(cscores_orig, cscores_new, tolerance = 1e-10)

if (isTRUE(equal_res)) {
    cat("[PASS] C-scores for ribo_toy are IDENTICAL.\n")
} else {
    cat("[FAIL] C-scores differ!\n")
    print(equal_res)
}

# Also verify MEAN method briefly
cat("\nVerifying 'mean' method on ribo_toy...\n")
cscores_orig_mean <- apply(counts, 2, compute_sample_cscore_orig, rna_factor = rnas, flanking = 6, method = "mean")
ribo_new_mean <- compute_cscore(ribo_toy, flanking = 6, method = "mean", ncores = 1)
cscores_new_mean <- assay(ribo_new_mean, "cscore")

equal_res_mean <- all.equal(cscores_orig_mean, cscores_new_mean, tolerance = 1e-10)
if (isTRUE(equal_res_mean)) {
    cat("[PASS] Mean C-scores for ribo_toy are IDENTICAL.\n")
} else {
    cat("[FAIL] Mean C-scores differ!\n")
    print(equal_res_mean)
}
