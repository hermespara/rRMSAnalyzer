# Helper script to update ribo_toy data
# Load dependencies
library(SummarizedExperiment)
library(S4Vectors)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)

# Source all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
sapply(r_files, source)

# Paths
count_path <- "inst/extdata/miniglioma"
metadata_path <- "inst/extdata/metadata.csv"

# Load data
ribo_toy <- load_ribodata(
  count_path = count_path,
  metadata = metadata_path,
  count_sep = "\t",
  metadata_sep = ",",
  count_header = FALSE,
  count_value = 3,
  count_rnaid = 1,
  count_pos = 2,
  metadata_key = "filename",
  metadata_id = "samplename",
  flanking = 6,
  method = "median",
  ncores = 1
)

# Function to rename RNA wrapper for ribo_toy if needed (usually done in examples)
# But standard ribo_toy should be raw loaded? 
# The original ribo_toy in tests seems to accept rename_rna in tests. 
# So raw load is fine.

save(ribo_toy, file = "data/ribo_toy.rda")
message("ribo_toy.rda updated.")
