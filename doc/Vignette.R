## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=F------------------------------------------
library(rXTalkViz)
library(dplyr)
library(data.table)

## ----message=FALSE, warning=F-------------------------------------------------
print(head(example_enrich.df))

## ----message=FALSE, warning=F-------------------------------------------------
filtered_data <- example_enrich.df %>%
  filter(p.adjust < 0.001)
filtered_data %>% nrow

## ----out.width="100%", out.height="100%", message=FALSE, warning=F------------
# filtered_df <- example_enrich.df %>%
#   filter(p.adjust < 0.001)
# xTalk_wrapper(filtered_df,
#               doPlot = T,
#               doXtalkQuant = T,
#               nPermute = 2,
#               min_cross_talk_score = 1.0,
#               plot_width = 10,
#               plot_height = 10)

## ----fig.width=10, fig.height=10, message=FALSE, warning=F--------------------
filtered_df <- example_enrich.df %>%
  filter(p.adjust < 0.001)

p <- xTalkPlot.Type_I.NetworkView(
  enrich.df = filtered_df,
  string_PPI_score_th = string_PPI_score_th,
  showCategory = 100,
  layout = "linear",
  colorEdge = TRUE,
  circular = TRUE)
p

