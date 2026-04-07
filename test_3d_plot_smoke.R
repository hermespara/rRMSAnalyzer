pkgload::load_all(".")
library(SummarizedExperiment)

# Load toy data
data("ribo_toy")

# Check if annotated
annots <- get_annotation(ribo_toy)
print("Initial annotations:")
print(annots)

# If no annotations, add some (ribo_toy usually has some if rename_rna was called)
# Let's ensure we have something to plot
if (nrow(annots) == 0) {
    print("Adding mock annotations...")
    # Mock annotation df
    mock_annot <- data.frame(
        rnapos = c(100, 200),
        rna = c("NR_046235.3_28S", "NR_046235.3_18S"),
        site = c("28S_TEST", "18S_TEST")
    )
    ribo_toy <- annotate_site(ribo_toy, mock_annot, anno_rna = "rna", anno_pos = "rnapos", anno_value = "site")
}

# Try to generate the plot object with performance optimizations
print("Testing plot_3d_structure function with simplification...")
p <- plot_3d_structure(ribo_toy, pdb_id = "4UG0", show_proteins = FALSE, cartoon_quality = 2)

print("Class of output:")
print(class(p))

# Check internal structure
print("Plot object structure (head elements):")
print(names(p$x))
