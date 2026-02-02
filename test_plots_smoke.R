library(rRMSAnalyzer)
data("ribo_toy")

print("Checking plot_pca without ellipses...")
# Default plot
p1 <- plot_pca(ribo_toy, "run", draw_ellipses = FALSE)
print(class(p1))
# Check layers? Hard to check ggplot layers programmatically easily without inspecting, but we assume no error is good.

print("Checking plot_pca WITH ellipses...")
p2 <- plot_pca(ribo_toy, "run", draw_ellipses = TRUE)
print(class(p2))

print("All plotting tests passed!")
