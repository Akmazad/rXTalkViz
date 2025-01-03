% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xTalk_wrapper.R
\name{xTalk_wrapper}
\alias{xTalk_wrapper}
\title{This function works as a wrapper for running the whole cross-talk analysis pipeline, including the Type-I and Type-II cross-talk quantification, visualization and the result exports.}
\usage{
xTalk_wrapper(
  enrich.df = NULL,
  showCategory = 100,
  doPlot = T,
  doXtalkQuant = T,
  layout = "linear",
  colorEdge = T,
  isCircular = T,
  nPermute = 10,
  neighbourhood_th = 1,
  string_PPI_score_th = 900,
  plot_width = 10,
  plot_height = 8,
  min_cross_talk_score = 0.1,
  outdir = "./xTalk_results/"
)
}
\arguments{
\item{enrich.df}{Enrichment object (a R dataframe) that has the following column names: "Description", "Pathway.geneSymbol","Overlapping.geneSymbol","Overlapping.geneID", "Count", "GeneRatio", "pvalue","p.adjust")}

\item{showCategory}{A number for thresholding the number pathways to be included in the analysis}

\item{doPlot}{A logical value to decide if the cross-talk visualization should be done while running the analysis}

\item{doXtalkQuant}{A logical value to decide if the cross-talk quantification should be done while running the analysis}

\item{layout}{A string that should indicate the layout name (see igraph layout names)}

\item{colorEdge}{A logical value to decide if the edges should be colored in the network plots}

\item{isCircular}{A logical value to decide if the network plot should be circular}

\item{nPermute}{A number to indicate how many times the permutations should be done for computing p-values for each Type-II cross-talks}

\item{neighbourhood_th}{A number to threshold the PPI neighborhood while calculating the network proximity}

\item{string_PPI_score_th}{A number to threshold for filtering PPI confidence scores from String Database}

\item{plot_width}{A number for specifying the plot width}

\item{plot_height}{A number for specifying the plot height}

\item{min_cross_talk_score}{A number to threshold the cross-talk scores to reduce the visulization}

\item{outdir}{A directory that indicates where the outputs will be saved (all the plot PDFs, raw ggplot2 objects and the final report as csv)}
}
\value{
""
}
\description{
This function works as a wrapper for running the whole cross-talk analysis pipeline, including the Type-I and Type-II cross-talk quantification, visualization and the result exports.
}
\examples{
library(rXTalkViz)
library(dplyr)
filtered_df <- example_enrich.df \%>\% filter(p.adjust < 0.001)
xTalk_wrapper(filtered_df,
              doPlot = TRUE,
              doXtalkQuant = TRUE,
              nPermute = 2,
              min_cross_talk_score = 1.0,
              plot_width = 10,
              plot_height = 10)
}
