#' A function to visualize (Network plot) Type-I cross-talk among pathways from enrichment analysis
#'
#' @param enrich.df Enrichment object (a R dataframe) that has the following column names: "Description", "Pathway.geneSymbol","Overlapping.geneSymbol","Overlapping.geneID", "Count", "GeneRatio", "pvalue","p.adjust")
#' @param showCategory A number for thresholding the number pathways to be included in the analysis
#' @param string_PPI_score_th A number to threshold for filtering PPI confidence scores from String Database
#' @param layout A string that should indicate the layout name (see igraph layout names)
#' @param colorEdge A logical value to decide if the edges should be colored in the network plots
#' @param circular A logical value to decide if the network plot should be circular
#'
#' @return  A ggplot2 object
#' @export ""
#'
#' @examples
#' library(rXTalkViz)
#' library(dplyr)
#' filtered_df <- example_enrich.df %>% filter(p.adjust < 0.001)
#' xTalkPlot.Type_I.NetworkView(enrich.df = filtered_df,
#'                             string_PPI_score_th = string_PPI_score_th,
#'                             showCategory = 100,
#'                             layout = "linear",
#'                             colorEdge = TRUE,
#'                             circular = TRUE)

xTalkPlot.Type_I.NetworkView <-
  function(enrich.df = NULL,
           string_PPI_score_th = 900,
           showCategory = 100,
           layout = "kk",
           colorEdge = T,
           circular = F) {

    if(is.null(enrich.df)){
      stop("Enrichment data has not been provided!")
    }
    # if(is.null(fc.dat)){
    #   warning("Fold-change data has not been provided!")
    # }

    n <- if(showCategory > nrow(enrich.df)) nrow(enrich.df) else showCategory
    geneSets <- str_split(enrich.df[, 3], "/") %>%
      `names<-`(enrich.df[,1]) %>%
      removeSoloNodes() # remove solo nodes that are not shared by minimum 2 pathways (unfit for being Type-I crosstalk)

    g <- geneSets %>%
      list2df()  %>%
      graph.data.frame(directed = F)


    # basically find the number of Enriched Pathways or Terms provided
    # this will set the node size: pathway node sizes are proportionate to their
    # constituent gene/protein number
    size <- sapply(geneSets, length)
    V(g)$size <- min(size) / 2
    n <- length(geneSets)
    n <- length(V(g)) - n
    V(g)$size[(n + 1):length(V(g))] <- size
    # color all the nodes
    V(g)$color <- rep("grey", length(V(g)))
    # now change the colors for the pathway nodes only
    V(g)$color[(n + 1):length(V(g))] <- rep("orange", length(geneSets))


    if (circular) {
      layout <- "linear"
      geom_edge <- geom_edge_arc
    }else{
      geom_edge <- geom_edge_link
    }

    # Modification on Type-II cross-talk ----
    # new_edges = data.frame()
    #
    # if(CrossTalk_Type == "Type_II"){
    #   for(i in 1:(length(geneSets)-1)){
    #     set_i = geneSets[i] %>% unlist(use.names = F)
    #     for(j in (i+1):length(geneSets)){
    #       set_j = geneSets[j] %>% unlist(use.names = F)
    #
    #       # extract PPI info between genes of different pathways only
    #       new_edges = protein.links %>%
    #         dplyr::filter(combined_score > string_PPI_score_th) %>%
    #         dplyr::select(-combined_score) %>%
    #         dplyr::filter((from %in% set_i & to %in% set_j) |
    #                         (from %in% set_j & to %in% set_i)) %>%
    #         dplyr::filter(!(from %in% set_i & to %in% set_i)) %>%
    #         dplyr::filter(!(from %in% set_j & to %in% set_j)) %>%
    #         rbind(new_edges, .)
    #
    #     }}
    #
    #   # remove duplicated rows, e.g., (i,j) and (j,i) are the same PPI
    #   new_edges <- new_edges[!duplicated(
    #     apply(new_edges,
    #           1,
    #           function(x) paste(sort(x), collapse = ''))),
    #     ]
    #   if(nrow(new_edges) > 0)
    #     g <- igraph::add_edges(g,
    #                            new_edges %>%
    #                              as.matrix() %>%
    #                              t %>%
    #                              c,
    #                            color = "red") %>% simplify()
    # }

    # -----
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    E(g)$alpha <- rep(0.15, length(E(g)))
    # edge_layer <- geom_edge(aes_(color = ~category, alpha = ~I(alpha)),
    #                         # color = "grey",
    #                         show.legend = F)
    edge_layer <- geom_edge(aes(alpha = stat(index)),
                            strength=0.2,
                            colour = 'darkgrey',
                            show.legend=F)

    # load Foldchange
    #
    #     if (!is.null(fc.dat)) { # fold-change data is provided
    #       colnames(fc.dat) = c("SYMBOL", "FC")
    #       temp = dplyr::inner_join(data.frame(SYMBOL = V(g)$name),
    #                                fc.dat, by = "SYMBOL")
    #
    #       palette <- fc_palette(-6:6)
    #       V(g)$color[1:n] <- temp$FC
    #       V(g)$degree = degree(g)
    #
    #     graph.obj <- ggraph(g,layout = layout, circular = circular) +
    #       edge_layer +
    #       geom_node_point(aes_(
    #         color =  ~ as.numeric(as.character(color)),
    #         # fill = "orange",
    #         size =  ~ size), alpha=0.2) +
    #       scale_color_gradientn(name = "Log2(Fold Change)",
    #                             colors = palette,
    #                             na.value = "#ffc34d") +
    #       scale_size(range = c(3, 10), breaks = unique(round(seq(
    #         min(size), max(size), length.out = 4)))) +
    #       theme_graph() +
    #       geom_node_text(aes_(label=~name), repel=TRUE)
    #     }else{
    # fold-change data is missing
    graph.obj <- ggraph(g,layout = layout, circular = circular) +
      edge_layer +
      geom_node_point(aes_(size =  ~ size,
                           fill =  ~ color), pch=21, show.legend = F) +
      geom_node_point(aes_(size =  ~ size,
                           fill =  ~ color,
                           alpha=0.2), pch=21, show.legend = F) +
      scale_size(range = c(3, 10), breaks = unique(round(seq(
        min(size), max(size), length.out = 4)))) +

      theme_graph() +
      geom_node_text(aes_(label=~name), repel=TRUE)
    # }

    return(graph.obj)

  }

#' A function to visualize (Sankey plot) Type-I cross-talk among pathways from enrichment analysis
#'
#' @param enrich.df Enrichment object (a R dataframe) that has the following column names: "Description", "Pathway.geneSymbol","Overlapping.geneSymbol","Overlapping.geneID", "Count", "GeneRatio", "pvalue","p.adjust")
#' @param outdir A directory that indicates where the outputs will be saved (all the plot PDFs, raw ggplot2 objects and the final report as csv)
#'
#' @return ""
#' @export ""
#'
#' @examples
#' ""
xTalkPlot.Type_I.Sankey <-function(enrich.df = NULL, outdir = NULL){

  if(is.null(enrich.df))
    stop("Enrichment Data can't be NULL !")
  if(is.null(outdir))
    stop("output directory is compulsory !")

  # library(networkD3)
  # require(stringi)
  # require(dplyr)
  # require(data.table)

  geneSets <- str_split(enrich.df[, 3], "/") %>%
    `names<-`(enrich.df[,1]) %>%
    removeSoloNodes() # remove solo nodes that are not shared by minimum 2 pathways (unfit for being Type-I crosstalk)


  pathway.data = names(geneSets)
  nTerms = length(pathway.data)

  enriched.genes = enrich.df[,3] %>%
    str_split(pattern = "/", simplify = F) %>%
    unlist(use.names = F) %>%
    unique()

  pathway.nodes = data.frame(name = pathway.data, group = "pathway")
  enriched.gene.nodes = data.frame(name = enriched.genes, group = "enriched.genes")
  nodes = bind_rows(pathway.nodes, enriched.gene.nodes)


  edges = data.frame(source = c(),
                     target = c(),
                     group = c())
  links = data.frame(source = c(),
                     target = c(),
                     group = c())

  for(i in 1:nTerms){
    source = geneSets[i] %>% names()
    target = geneSets[[i]]
    group = "enriched.genes"

    edges = rbind(edges, cbind(source, target, group))
  }
  nodes = inner_join(nodes, data.frame(
    name = c(edges[, 1] %>% as.character(),
             edges[, 2]  %>% as.character()) %>% unique()), by = "name")


  links = data.frame(
    source = base::match(edges[, 1] %>% as.character(), nodes[, 1])-1,
    target = base::match(edges[, 2] %>% as.character(), nodes[, 1])-1,
    # group = base::match(edges[, 3] %>% as.character(), nodes[, 2])-1
    value = 1
  )

  p <- sankeyNetwork(Links = links, Nodes = nodes, Source = 'source', Value = 'value', fontFamily = "Helvetica",
                     Target = 'target', NodeID = 'name', sinksRight = T, nodePadding = 0,
                     fontSize = 8, nodeWidth = 10, iterations = 0)

  # library(htmlwidgets)
  saveWidget(p, file=paste0(outdir,"xTalkSankey_Type_I.html"))
}

# xTalkPlot.Netmodel <-
#   function(enrich.df,
#            expr.dat = Null,
#            showCategory = 100,
#            fc.dat   = NULL,
#            layout = "kk",
#            colorEdge = T,
#            circular = F,
#            node_label = "all"){
#     if(is.null(expr.dat)){
#       stop("Expression dat can't be null")
#     }
#     if(class(expr.dat) != "data.frame") expr.dat = expr.dat %>% as.data.frame()
#
#   }


#' A function to visualize (Network plot) Type-II cross-talk among pathways from enrichment analysis
#'
#' @param df A R dataframe that has pair-wise cross-talk information (measured via network proximity scores)
#' @param isCircular A logical value to decide if the network plot should be circular
#' @param layout A string that should indicate the layout name (see igraph layout names)
#'
#' @return A ggplot2 object
#' @export ""
#'
#' @examples ""
#' ""

xTalkPlot.Type_II.NetworkView <- function(df = NULL, isCircular, layout){

  if(is.null(df)){
    stop("df can't be null")
  }
  minEdgeWidth = min(df$weight, na.rm = T)
  maxEdgeWidth = max(df$weight, na.rm = T)


  g <- graph_from_data_frame(df)
  V(g)$degree = degree(g)
  V(g)$frame.color <- "white"
  V(g)$color <- "orange"
  E(g)$weight <- df$weight

  # E(g)$category <- rep(names(V(g)), sapply(V(g), length))
  # print(E(g)$category)
  E(g)$alpha <- rep(0.15, length(E(g)))

  if (isCircular) {
    layout <- "linear"
    geom_edge <- geom_edge_arc
  }else{
    geom_edge <- geom_edge_link
  }


  edge_layer <- geom_edge(aes_(
    # width = ~ weight,
    alpha = ~I(alpha)
    # alpha = after_stat(index)
  ),
  show.legend = F)

  graph.obj <- g %>% ggraph(
    layout = layout,
    circular = isCircular) +
    edge_layer +
    geom_node_point(aes_(size =  ~ degree,
                         fill =  ~ color), pch=21, show.legend = F) +
    # geom_edge_arc(aes(alpha = after_stat(index)), strength = 0.2) +
    geom_node_point(aes_(size =  ~ degree,
                         fill =  ~ color,
                         alpha=0.2), pch=21, show.legend = F) +
    # scale_edge_width(range = c(minEdgeWidth, maxEdgeWidth))+ # control size
    theme_graph() +
    theme(legend.position = "none") +
    geom_node_text(aes_(label=~name), repel=TRUE)


  return(graph.obj)

}


#' A function to visualize (Heatmap plot) Type-II cross-talk among pathways from enrichment analysis
#'
#' @param df A R dataframe that has pair-wise cross-talk information (measured via network proximity scores)
#'
#' @return a heatmap object from ComplexHeatmap package
#' @export ""
#'
#' @examples ""
#' ""

xTalkPlot.Type_II.heatmap <- function(df = NULL){

  if(is.null(df)){
    stop("df can't be null")
  }


  plot.df <- df %>%
    pivot_wider(names_from = Pathway_2,
                values_from = `weight`) %>%
    column_to_rownames("Pathway_1")
  plot.df[plot.df == Inf] <- 0
  plot.df <- plot.df %>%
    na.replace(fill = 0) %>%
    as.matrix()

  # browser()
  col_fun = colorRamp2(c(min(plot.df, na.rm = T),
                         max(plot.df, na.rm = T)),
                       c("white", "red"))

  # browser()
  lgd = Legend(col_fun = col_fun,
               title = "Type-II cross-talk",
               direction = "horizontal",
               title_position  = "topcenter")

  ht <- plot.df %>%
    Heatmap(
      # heatmap_legend_param = list(
      #   title = "Type-II cross-talk",
      #   direction = "horizontal",
      #   title_position  = "topcenter"
      #   ),
      # col = col_fun
    )

  draw(ht,
       show_heatmap_legend = F
       # heatmap_legend_side="bottom"
  )
  # draw(lgd, x = unit(0.95, "npc"), y = unit(1, "cm"), just=c("right","bottom"))

}

update_n <- function(df, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  n <- showCategory
  if (nrow(df) < n) {
    n <- nrow(df)
  }

  return(n)
}

extract_geneSets <- function(df, n) {
  n <- update_n(df, n)
  df <- df[1:n, ]
  geneSets <- str_split(df[, 3], "/")
  names(geneSets) = df[, 1]
  return(geneSets) ## if n is a vector of Description
}

list2graph <- function(inputList) {
  x <- list2df(inputList)
  # print(x)
  g <- graph.data.frame(x, directed = FALSE)
  return(g)
}

fc_palette <- function(fc) {
  if (all(fc > 0, na.rm = TRUE)) {
    palette <- enrichplot::color_palette(c("gray", "red"))
  } else if (all(fc < 0, na.rm = TRUE)) {
    palette <- enrichplot::color_palette(c("darkgreen", "gray"))
  } else {
    palette <-
      enrichplot::color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
    # enrichplot::color_palette(c("darkgreen", "gray", "red"))

  }
  return(palette)
}

removeSoloNodes <- function(aList){
  df <- aList %>% list2df() %>%
    `colnames<-`(c("source", "target"))

  # retrieve the genes that ware only enriched by one pathway, hence
  # unfit to be defined as Type-I cross-talk
  soloGenes <- df %>%
    dplyr::count(source) %>%
    dplyr::filter(n == 1) %>%
    dplyr::select(source) %>%
    unlist(use.names = F)

  # exclude the sologenes
  df <- df %>%
    dplyr::filter(! source %in% soloGenes)
  tapply(df$source, df$target, list) %>%
    return()
}
