# Inspect ribo_toy
library(SummarizedExperiment)

load("data/ribo_toy.rda")
print(class(ribo_toy))
print(colnames(colData(ribo_toy)))
print(colData(ribo_toy))
