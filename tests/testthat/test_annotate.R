test_that("Annotation works correctly", {
  data("ribo_toy")

  # We will create a custom annotation just for this test

  custom_anno <- data.frame(
    rnapos = c(15, 76, 100, 401),
    rna = c("5.8S", "5.8S", "18S", "28S"),
    nomenclature = c("5.8S_Um14", "5.8S_Gm75", "18S_Am99", "28S_Am391")
  )

  ribo_toy <- rename_rna(ribo_toy)
  ribo_toy <- annotate_site(ribo_toy, custom_anno)
  ribom <- extract_data(ribo_toy, only_annotated = T)

  expected_ribom <- readRDS(testthat::test_path("testdata", "anno.rds"))

  # Align rows by site name for robust comparison
  ribom <- ribom[match(expected_ribom$site, ribom$site), ]

  testthat::expect_equal(ribom, expected_ribom, ignore_attr = TRUE)
})
