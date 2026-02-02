# Helper script to run tests locally without full package installation
# Load dependencies
library(SummarizedExperiment)
library(S4Vectors)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(testthat)
library(cli)
# Add other dependencies if needed

# Source all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
sapply(r_files, source)

# Run tests
testthat::test_dir("tests/testthat")
