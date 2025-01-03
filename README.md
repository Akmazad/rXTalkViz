# rXTalkViz: A tool for quantification and visualization of functional cross-talks from enrichment analysis
Original article can be cited as:
```
@article{Azad2024,
  title = {rXTalkViz: A R package to quantify,  visualize,  and report carcinogenic footprints of functional pathway cross-talks},
  volume = {22},
  ISSN = {2665-9638},
  url = {http://dx.doi.org/10.1016/j.simpa.2024.100708},
  DOI = {10.1016/j.simpa.2024.100708},
  journal = {Software Impacts},
  publisher = {Elsevier BV},
  author = {Azad,  A.K.M.},
  year = {2024},
  month = nov,
  pages = {100708}
}
```
## Introduction
The `rXTalkViz` package provides tools for quantifying and visualizing cross-talks between biological pathways. It models the pathway cross-talks in two categories:

1. **Type-I cross-talks**: Via shared nodes (genes/proteins/molecules) - quantified with Jaccard-Index score.
2. **Type-II cross-talks**: Via shared interactions (Protein-Protein interactions, or data-driven interactions) - quantified with Network proximity scores.

Visualizations include network plots, Sankey diagrams, and heatmap plots. The quantification results are also saved as a spreadsheet (CSV file), which can be an input for other bioinformatics pipeline. The plot objects (ggplot2) are also saved within the result directory so that the users can visualize them with their own preferred settings. This guide will walk you through the basic usage of the package, including data loading, cross-talk analysis, and result visualization.

## Installation
```r
library(devtools)
install_github("Akmazad/rXTalkViz")
```

### Step 1: Loading the Required Libraries
First, load the necessary libraries, including `rXTalkViz` and other dependencies such as `dplyr` and `data.table`.

```r
library(rXTalkViz)
library(dplyr)
library(data.table)
```
### Step 2: Loading Example Data
The package comes with an example dataset named example_enrich.df. Let's observe the head of this data to understand its structure.

```r
print(head(example_enrich.df))
```

This dataset contains enrichment results for various pathways. We will use this data to explore the cross-talk analysis features of the package.

### Step 3: Filtering Data
To focus on the most significant pathways, we will filter the dataset to include only those pathways with an adjusted p-value (p.adjust) below 0.001. Moreover, this will greatly control the time required for running the whole pipeline by invoking the wrapper function.


```r
filtered_df <- example_enrich.df %>%
  filter(p.adjust < 0.001)
filtered_df %>% nrow
```

### Step 4: Performing Cross-Talk Analysis
With the filtered data, we can now perform cross-talk analysis using the `xTalk_wrapper` function. This function allows you to quantify, visualize cross-talks between pathways, and save all the results within the specified directory. Note, if the directory doesn't exist, it will create on but with a warning, and the directory name can be specified by changing the 'outdir' parameter.

```r
xTalk_wrapper(filtered_df,
             doPlot = TRUE,
             doXtalkQuant = TRUE,
             nPermute = 2,
             min_cross_talk_score = 1.0,
             plot_width = 10,
             plot_height = 10)
```

```r
filtered_df <- example_enrich.df %>%
  filter(p.adjust < 0.001)

p <- xTalkPlot.Type_I.NetworkView(
  enrich.df = filtered_df,
  string_PPI_score_th = string_PPI_score_th,
  showCategory = 100,
  layout = "linear",
  colorEdge = T,
  circular = T)
p
```

In this network plot, Type-I cross-talks is shown among highly enriched pathways via shared genes.

## Conclusion
In this vignette, we demonstrated the basic functionality of the rXTalkViz package. We loaded example data, filtered it based on significance, and performed cross-talk analysis. For more detailed analysis, refer to the package documentation.

## Acknowledgement
We thank Dr. Guangchuang (YuLabSMU) for his excellent R package: clusterprofiler, as some parts of the network visualization code in our work have been reproduced/inspired by that package. 
