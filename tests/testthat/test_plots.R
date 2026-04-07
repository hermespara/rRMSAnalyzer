test_that("test if plots creation do not fail", {
  data(ribo_toy)

  custom_anno <- data.frame(
    rnapos = c(15, 76, 100, 401),
    rna = c("5.8S", "5.8S", "18S", "28S"),
    nomenclature = c("5.8S_Um14", "5.8S_Gm75", "18S_Am99", "28S_Am391")
  )

  ribo_toy <- rename_rna(ribo_toy)
  ribo_toy <- annotate_site(ribo_toy, custom_anno)

  expect_no_error(plot_pca(ribo_toy, "condition"))
  expect_no_error(plot_pca(ribo_toy, "condition",
    sites = c("5.8S_Um14", "18S_Am99")
  ))
  expect_error(
    plot_pca(ribo_toy, "condition", sites = c("FAKE_SITE")),
    "No matching sites"
  )
  expect_error(
    plot_pca(ribo_toy, "condition",
      only_annotated = TRUE,
      sites = c("5.8S_Um14")
    ),
    "cannot be used together"
  )
  expect_warning(
    plot_pca(ribo_toy, "condition", sites = c("5.8S_Um14", "18S_Am99", "FAKE")),
    "not found"
  )
  expect_no_error(plot_coa(ribo_toy, "condition"))
  expect_s3_class(
    plot_boxplot_count(
      ribo_toy,
      color_col = "condition",
      sample_colors = c("RNA ref" = "grey40", "cond1" = "blue", "cond2" = "red")
    ),
    "ggplot"
  )
  expect_error(
    plot_boxplot_count(
      ribo_toy,
      color_col = "condition",
      sample_colors = c("RNA ref" = "grey40")
    ),
    "Missing color mapping"
  )

  expect_no_error(plot_boxplot_count(ribo_toy))
  expect_no_error(plot_boxplot_cscores(ribo_toy))
  expect_no_error(plot_rle(ribo_toy))

  expect_s4_class(plot_heatmap(ribo_toy,
    color_col = c("run", "condition")
  ), "Heatmap")
  expect_s4_class(plot_heatmap_corr(ribo_toy, "count", "run"), "Heatmap")

  expect_no_error(plot_counts_env(ribo_toy, "5S", 50))
  expect_no_error(plot_counts_env(ribo_toy, "5S", 50, c("S1", "S2")))

  expect_no_error(plot_diff_sites(ribo_toy, "condition"))
})
