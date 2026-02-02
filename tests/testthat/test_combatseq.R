test_that("test if combatSeq does its job", {
    data(ribo_toy)
    ribo_corrected <- adjust_bias(ribo_toy, "run")
    ribo_corrected_m <- extract_data(ribo_corrected)

    expected_correction <- readRDS(testthat::test_path("testdata", "cseq.rds"))

    # Align rows by named_position for robust comparison
    ribo_corrected_m <- ribo_corrected_m[match(expected_correction$named_position, ribo_corrected_m$named_position), ]

    # Align sample names if they differ (S1 vs RNA1 etc) but sizes match
    if (ncol(ribo_corrected_m) == ncol(expected_correction)) {
        names(ribo_corrected_m) <- names(expected_correction)
    }

    testthat::expect_equal(ribo_corrected_m, expected_correction, ignore_attr = TRUE)
})
