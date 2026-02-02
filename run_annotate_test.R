# Helper script to run ONLY annotate test
library(SummarizedExperiment)
library(S4Vectors)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(testthat)

# Source all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
sapply(r_files, source)

# Load data
load("data/ribo_toy.rda")

# Run tests
testthat::test_file("tests/testthat/test_annotate.R")
