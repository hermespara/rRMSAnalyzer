library(renv)
renv::init()

renv::install("hermespara/rRMSAnalyzer")

# installation du package en 15 minutes

library(rRMSAnalyzer)

data("ribo_toy")
data("human_methylated")
ribo_toy <- rename_rna(ribo_toy)
ribo_toy <- annotate_site(ribo_toy,human_methylated)
plot_sites_by_IQR(ribo = ribo_toy, plot = "IQR")

ribo_toy <- compute_cscore(ribo = ribo_toy)

plot_pca(ribo = ribo_toy, color_col = "condition", only_annotated = TRUE)

plot_coa(ribo = ribo_toy, color_col = "condition", only_annotated = TRUE)

count_data <- extract_data(ribo = ribo_toy, col = "counts") 

plot_counts_env(ribo = ribo_toy, rna = "5.8S", pos = 14)

plot_counts_fraction(ribo = ribo_toy)

plot_diff_sites(ribo = ribo_toy, factor_column = "condition", p_cutoff = 0.1)
