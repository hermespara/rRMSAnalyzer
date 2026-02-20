<!-- README.md is generated from README.Rmd. Please edit that file -->

# rRMSAnalyzer: Comprehensive Analysis of rRNA 2'-O-Methylation

**rRMSAnalyzer** is an R package designed for the analysis and visualization of **rRNA 2'-O-ribose methylation (2'-O-Me)**, a critical chemical modification that fine-tunes ribosomal function and translation accuracy.

Using data from **RiboMethSeq**, an RNA-seq-based protocol, rRMSAnalyzer enables researchers to quantify 2'-O-Me levels with high precision. By computing **C-scores** (methylation scores) from read-end counts, the package provides a robust framework for exploring methylation landscapes across different biological conditions.

### Key Features (Version 2)

*   **Precise Quantification**: Compute C-scores using flexible parameters (mean or median local coverage).
*   **Batch Effect Correction**: Integrated **ComBat-Seq** support to remove technical biases and ensure robust comparisons.
*   **Rich Visualization**: Generate publication-ready plots (PCA, Heatmaps, Boxplots) to explore methylation patterns.
*   **Annotated Reference Data**: Includes curated lists of known human rRNA methylation sites.
*   **Automated Reporting**: Generate comprehensive quality control reports with a single function.
*   **Seamless Integration**: Built on the `SummarizedExperiment` class for interoperability with the Bioconductor ecosystem.

> **Note**: For processing raw sequencing data (FASTQ) into read-end counts compatible with rRMSAnalyzer, we recommend our dedicated [Nextflow pipeline](https://github.com/RibosomeCRCL/ribomethseq-nf).

## Installation

The latest version of rRMSAnalyzer package can be installed from Github
with:

``` r
library(devtools)
devtools::install_github("RibosomeCRCL/rRMSAnalyzer")
```

## Usage

``` r
library(rRMSAnalyzer)

ribo <- load_ribodata(
              count_path = "/path/to/your/csvfiles/directory/",
              metadata = "path/to/metadata.csv",
              metadata_key = "filename",
              metadata_id = "samplename")

# Compute the c-score using different parameters,
# including calculation of the local coverage using the mean instead of the median
ribo <- compute_cscore(ribo, method = "mean")

# If necessary, adjust any technical biases using ComBat-Seq.
# Here, as an example, we use the "library" column in metadata.
ribo <- adjust_bias(ribo,"library")

# Plot a Principal Component Analysis (PCA) whose colors depend on the "condition" column in metadata
plot_pca(ribo,"condition")
```

## Documentation

For a comprehensive guide on using **rRMSAnalyzer**, please visit our [official website](https://ribosomecrcl.github.io/rRMSAnalyzer/).

*   [Getting Started](https://ribosomecrcl.github.io/rRMSAnalyzer/articles/rRMSAnalyzer.html): A step-by-step tutorial.
*   [Reference](https://ribosomecrcl.github.io/rRMSAnalyzer/reference/index.html): Detailed function documentation.

*A toy dataset (`ribo_toy`) is included in the package for testing and demonstration.*

## Support & Contributions

We welcome feedback and contributions!
*   **Report Bugs or Suggest Features**: Please open an issue on our [GitHub Issues page](https://github.com/RibosomeCRCL/rRMSAnalyzer/issues).
*   **Contribute**: Feel free to submit puly requests to improve the package.

## Acknowledgements

We would like to thank all our collaborators from Jean-Jacques Diaz Team
and the Bioinformatic Platform Gilles Thomas for their advices and
suggestions.

## Funding

This project has been funded by the French Cancer Institute (INCa, PLBIO
2019-138 MARACAS), the SIRIC Program (INCa-DGOS-Inserm_12563 LyRICAN),
LabEX program (DEVweCan), the French association Ligue Nationale Contre
le Cancer and Synergie Lyon Cancer Foundation.
