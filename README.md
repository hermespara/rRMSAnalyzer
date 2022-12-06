<!-- README.md is generated from README.Rmd. Please edit that file -->

# Installation

You can install the latest rRMSAnalyzer version on CRAN with :

``` r
install.packages("rRMSAnalyzer")
```

# Usage

``` r
library(rRMSAnalyzer)

ribo <- load_ribodata("/path/to/your/csvfiles/directory/",
              "path/to/metadata.csv",
              metadata_key = "filename",
              metadata_id = "samplename")

#recompute c-score using window's mean instead of median.
ribo <- compute_cscore(ribo, method = "mean")

#Correct any technical biases with ComBat-Seq, using the library column in metadata
ribo <- adjust_bias(ribo,"library")

# Make a PCA plot and color it by biological condition (column in metadata)
plot_PCA(ribo,"biological_condition")
```

    #> [SUCCESS] Your data have been imported and the following RiboClass has been created :
    #> a RiboClass with 14 samples and 4 RNA(s) :
    #> Name : NR_023363.1_5S, length : 121 
    #> Name : NR_046235.3_5.8S, length : 157 
    #> Name : NR_046235.3_18S, length : 1869 
    #> Name : NR_046235.3_28S, length : 5070

# Help, bug reports and suggestions

If you find a bug, or have a suggestion to improve the package, open a new issue on : *github issue link coming soon!*

# Acknowledgement

We would like to thank

# Funding

This project has been funded.

# Input data

## The RiboClass

The RiboClass is the object used by this package to contain both your data and metadata. It is created when calling load_ribodata (see [#quickstart](#quickstart) or loading data).

It contains three main parts :

1.  **Data** : a list of dataframe, containing your samples’ count + calculated c-score.

2.  **Metadata** : a dataframe containing all informations related to your samples.

3.  **RNA_names** : current and original names of your RNAs in your data.

Other elements in the riboClass serve to keep a tab on some function calls’ parameters.

## GenomeCov-like Count data

To use this package, you need the following columns :

1.  The name of the RNA on which the 5’/3’ end count has been done on

2.  The position on the RNA

3.  The 5’/3’ count on this position

You can see a small example below :

| RNA | Position on RNA | 5’end count |
|-----|-----------------|-------------|
| 18S | 123             | 3746        |
| 18S | 124             | 345         |
| 18S | 125             | 324         |
| 18S | 126             | 789         |
| 18S | 127             | 1234        |

(it is not necessary to have an header in your count files. You can give the column numbering if there is no header).

The folder structure can be whatever you wish, as long as either the directory and its sub-directories contain the necessary CSV files.

**If you do not specify metadata, rRMSAnalyzer will try to get any CSV files.**

## Metadata

Only two columns are necessary for your metadata :

1.  filename : name of the csv file on disk. **Do not modify it unless the filename has changed on disk.**

2.  samplename : sample name you want to see on your plots and analyses. You are free to modify this column, as long as your names are unique.

You can then add any column in your metadata.

If you do not give any metadata in load_ribodata, an empty one will be created with filename and samplename columns pre-completed.

## Import new sites annotation

When you first load your data into a RiboClass, your sites are not annotated. But maybe you want to annotate some of them. You can use the [Annotations included with this package](#annotations-included-with-this-package), or create your own :

You need a dataframe with the following informations :

1.  RNA name, **matching the RNA name in your RiboClass.**

2.  Position on RNA.

3.  Name you want to give to the site.

You can see an example below :

| Position | rRNA | Nomenclature |
|----------|------|--------------|
| 15       | 5.8S | Um14         |
| 76       | 5.8S | Gm75         |
| 28       | 18S  | Am27         |

Example of annotation in human_methylated

Other columns in your annotation table will be ignored.

You can annotate you sites with [Annotate sites](#annotate-sites).

## Annotations included with this package

This package includes two dataframes with sites annotation :

-   human_methylated, containing the 112 known methylation sites for Human.

-   human_suspected, containing the 17 sites that are possible methylation sites.

You can annotate you sites with [Annotate sites](#annotate-sites).

# C-score calculation

## What is C-score

This is a score that has been introduced by Birkedal et al. It corresponds to the 2’Ome level at a rRNA position known to be methylated. The C-score represents a drop in end read coverage at a given position compared to the environmental coverage as described by Birkedal et al, 2015. This score can have a value between **0** (never methylated) and **1** (always methylated).

Two different C-score can be found in literature, depending if the authors the median or the mean of the window’s values. Both are implemented in this package.

## C-score computation when loading data

When load_ribodata is used, a c-score is automatically computed for all positions.

You can modify directly the parameters related to this computation:

``` r
load_ribodata("/path/to/csv/",
              "/path/to/metadata.csv",
              # everything below is linked to c-score computation
              flanking = 6, # window size
              method = "median", # use mean or median on window's values
              ncores = 8 # number of CPU cores to use for computation
              )
```

## Update C-Score

If you want to modify the how the c-score is computed and update the c-score in a RiboClass, you can use compute_cscore function.

In the following example, we modify both the window’s size and the computation method (mean instead of median) :

``` r
ribo <- compute_cscore(ribo,
                       flanking = 8,
                       method = "mean")
```

It will return a RiboClass with the new c-score. Please note that it will override the previous c-score.

# Quality control

## Samples PCA for all positions

We can spot a technical bias when we make a PCA of c-score for all genomic positions and color it by library (or batch).

(see also [Visualization with PCA](#visualization-with-pca) for more usages)

``` r
plot_PCA(ribo,"run")
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%"/>

The sample should not be grouped by library or batch. Unfortunately, this is clearly visible on PC 1. The following section will resolve this issue.

## Batch effect adjustment

If the sequencing has been done in several batches, a batch effect can occur. One way to remove this effect is by using ComBat-Seq adjustment. The rRMSAnalyzer package has a wrapper to use this adjustment on a RiboClass : *adjust_bias*.

``` r
ribo_adjusted <- adjust_bias(ribo,"run")
#> Found 2 batches
#> Using null model in ComBat-seq.
#> Adjusting for 0 covariate(s) or covariate level(s)
#> Estimating dispersions
#> Fitting the GLM model
#> Shrinkage off - using GLM estimates for parameters
#> Adjusting the data
#> Recalculating c-score using previously given parameters...
#> c-score method : median 
#> flanking window : 6
```

This function will return a new RiboClass with adjusted count values. The c-score will be automatically recomputed using the previous parameters given.

Now, if we redo the PCA, this time on our ComBat-Seq-adjusted RiboClass :

``` r
plot_PCA(ribo_adjusted,"run")
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%"/>

The bias on PC 1 is no longer visible !

## Keep or remove samples

If you want to analyze only a subset of samples, you can create a new RiboClass with keep_ribo_samples. In this example, we keep only the samples “2400” and “3307” :

``` r
ribo_2samples <- keep_ribo_samples(ribo_adjusted,c("938","2633"))
print(ribo_2samples)
#> a RiboClass with 2 samples and 4 RNA(s) :
#> Name : NR_023363.1_5S, length : 121 
#> Name : NR_046235.3_5.8S, length : 157 
#> Name : NR_046235.3_18S, length : 1869 
#> Name : NR_046235.3_28S, length : 5070
```

On the contrary, you can create a new RiboClass by removing some samples. Here, we remove “2400” and “3307” :

``` r
ribo_removed_samples <- remove_ribo_samples(ribo,c("938","2633"))
print(ribo_removed_samples)
#> a RiboClass with 12 samples and 4 RNA(s) :
#> Name : NR_023363.1_5S, length : 121 
#> Name : NR_046235.3_5.8S, length : 157 
#> Name : NR_046235.3_18S, length : 1869 
#> Name : NR_046235.3_28S, length : 5070
```

in both cases, only the remaining samples’ metadata are kept in the RiboClass object. You do not need to update it by hand.

# Data exploration

## Remove RNA

You can remove a RNA among all your RiboClass’ samples if it is not necessary for your analysis. Here, we do not need the 5S rRNA.

``` r
ribo_adjusted <- remove_rna(ribo,"NR_023363.1_5S")
print(ribo_adjusted)
#> a RiboClass with 14 samples and 3 RNA(s) :
#> Name : NR_046235.3_5.8S, length : 157 
#> Name : NR_046235.3_18S, length : 1869 
#> Name : NR_046235.3_28S, length : 5070
```

## Rename RNA

For some plots or functions, the RNA names can be important. **RNA names must match for sites annotation.**

Before annotating our sites using human_methylated, let’s check if our RiboClass’s RNA names are the same as those in this dataset :

``` r
cat("human_methylated's rna names :",unique(human_methylated$rRNA),"\n")
#> human_methylated's rna names : 5.8S 18S 28S
cat("ribo's rna names :",as.character(ribo_adjusted$rna_names$current_name))
#> ribo's rna names : NR_046235.3_5.8S NR_046235.3_18S NR_046235.3_28S
```

Our names are different ! We need to update them before annotation.

**We can use rename_rna**. We give the new RNA names by size order :

``` r
ribo_adjusted <- rename_rna(ribo_adjusted, c("5.8S","18S","28S"))
```

## Annotate RNA sites

Among the thousands sites belonging to the rRNA, some of them are more particular than others. In our case, we are interested by 2’Ome sites in the human model.

We can **annotate these sites with the included *human_methylated* dataset** :

``` r
ribo_adjusted <- annotate_site(ribo_adjusted,
                               annot = human_methylated,
                               anno_rna = "rRNA",
                               anno_pos = "Position",
                               anno_value = "Nomenclature")
```

If there is a problem with the RNA names, you should check [Rename RNA](#rename-rna).

This vignette has also some explanations on how to create your own sites annotation dataset with [Import new sites annotation](#import-new-sites-annotation).

## Visualization with PCA

After annotating these sites, it would be interesting to plot them. This is possible with plot_PCA function, using *only_annotated* parameter. For our example, we will also color it by biological condition :

``` r
plot_PCA(ribo_adjusted,
         color_col = "condition",
         only_annotated = T)
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%"/>

By default, the axes for your PCA plot are PC 1 and PC 2. **The *axes* parameter lets you choose which principal components you want to plot.** Using our previous example, we will plot this time PC 1 and PC 3 :

``` r
plot_PCA(ribo_adjusted,
         color_col = "condition",
         axes = c(1,3),
         only_annotated = T)
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="100%"/>

If you to do more than what is given by plot_PCA, **you can return the full dudi.pca object instead of a ggplot** by setting *pca_object_only* to True :

``` r
pca <- plot_PCA(ribo_adjusted,
         color_col = "condition",
         only_annotated = T,
         pca_object_only = T)
```

## Export data

If you want to use your results outside of this package, there are two ways to export your data :

### Export as a single matrix

You can use extract_data to export your data as a single dataframe

``` r
ribo_df <- extract_data(ribo_adjusted,
                        col = "cscore")
```

This will export all positions. If you want to export only annotated positions, you can set *only_annotated* parameter to True.

``` r
ribo_df <- extract_data(ribo_adjusted,
                        col = "cscore",
                        only_annotated = T)
```

### Export as a ggplot-ready dataframe

**If you want a ggplot-ready dataframe**, you can use format_to_plot. With the RiboClass, you can add additional metadata column if necessary.

``` r
ggplot_table <- format_to_plot(ribo_adjusted,"condition")
knitr::kable(ggplot_table[501:510,],caption = "excerpt from the output dataframe")
```

|     | sampleID | Cscore            | condition |
|:----|:---------|:------------------|:----------|
| 501 | 2346     | 0.946129580198981 | cond2     |
| 502 | 2346     | 0                 | cond2     |
| 503 | 2346     | 0                 | cond2     |
| 504 | 2346     | 0.203622392974753 | cond2     |
| 505 | 2346     | 0.295647558386412 | cond2     |
| 506 | 2346     | 0.552547770700637 | cond2     |
| 507 | 2346     | 0.887132352941176 | cond2     |
| 508 | 2346     | 0.810690423162584 | cond2     |
| 509 | 2346     | 0                 | cond2     |
| 510 | 2346     | 0.349796548337309 | cond2     |

excerpt from the output dataframe

### Export as a dataframe by condition

This export can be used to compare the mean count or c-score for each position between conditions (or other metadata).

``` r
mean_tb <- mean_samples_by_conditon(ribo_adjusted,"cscore","condition")
knitr::kable(mean_tb[501:510,],caption = "excerpt from the output dataframe")
```

| named_position | condition |      mean |        sd |
|:---------------|:----------|----------:|----------:|
| 18S_0167       | RNA ref   | 0.9826490 | 0.0033169 |
| 18S_0168       | cond1     | 0.1453841 | 0.0840809 |
| 18S_0168       | cond2     | 0.1762492 | 0.0405325 |
| 18S_0168       | RNA ref   | 0.0244490 | 0.0345761 |
| 18S_0169       | cond1     | 0.0000000 | 0.0000000 |
| 18S_0169       | cond2     | 0.0000000 | 0.0000000 |
| 18S_0169       | RNA ref   | 0.0000000 | 0.0000000 |
| 18S_0170       | cond1     | 0.0000000 | 0.0000000 |
| 18S_0170       | cond2     | 0.0000000 | 0.0000000 |
| 18S_0170       | RNA ref   | 0.0000000 | 0.0000000 |

excerpt from the output dataframe

# Reference
