#' A function to quantify Type-I cross-talk among pathways from enrichment analysis
#'
#' @param neighbourhood_th A number to threshold the PPI neighborhood while calculating the network proximity
#' @param string_PPI_score_th A number to threshold for filtering PPI confidence scores from String Database
#' @param pathway_a A set of strings which indicates gene symbols
#' @param pathway_b A set of strings which indicates gene symbols
#'
#' @return a vector of values: 1) the jaccard-index (as Type-I cross-talk score), 2) number of common/shared genes, 3) the common gene symbols (comma seperated)
#' @export ""
#'
#' @examples ""
#' ""
quantify_crossTalk_jaccard_index <- function(
    neighbourhood_th,
    string_PPI_score_th,
    pathway_a = c(),
    pathway_b = c()
    # pathway_a = c("EGFR", "FGFR3", "PIK3CA","AKT3", "MDM2","TP53"),
    # pathway_b=c("EGFR", "FGFR3", "PIK3CA","AKT3", "MDM2","XX")
){
  # see above note
  nCommons = intersect(pathway_a, pathway_b) %>% length()
  nAll = union(pathway_a, pathway_b) %>% length()

  ji_score = nCommons/nAll
  return(c(ji_score,
           nCommons,
           intersect(pathway_a, pathway_b) %>%
             paste0(collapse = ";")))
}
# quantify_crossTalk_MLE <- function(
#     neighbourhood_th,
#     pathway_a=c(),
#     pathway_b=c()){
#   # DONATO ET AL: https://genome.cshlp.org/content/23/11/1885
#
#   # load stringdb
#   string_db <- STRINGdb$new( version="11.5", species=9606,
#                              score_threshold=200, network_type="full", input_directory="")
# }

# quantify_crossTalk_functional_proximity <- function(
#     neighbourhood_th,
#     string_PPI_score_th,
#     pathway_a = c(),
#     pathway_b = c(),
#     nPermute
# ){
#
# }

#' A function to quantify Type-II cross-talk among pathways from enrichment analysis
#'
#' @param neighbourhood_th A number to threshold the PPI neighborhood while calculating the network proximity
#' @param string_PPI_score_th A number to threshold for filtering PPI confidence scores from String Database
#' @param pathway_a A set of strings which indicates gene symbols
#' @param pathway_b A set of strings which indicates gene symbols
#' @param nPermute A number to indicate how many times the permutations should be done for computing p-values for each Type-II cross-talks
#'
#' @return Network proximity score as the Type-II cross-talk value
#' @export ""
#'
#' @examples ""
#' ""
quantify_crossTalk_network_proximity <- function(
    neighbourhood_th,
    string_PPI_score_th,
    pathway_a = c(),
    pathway_b = c(),
    nPermute
    ){

  data(protein.links)
  # %>% as_tibble()
  #   dplyr::mutate(c("combined_score"), as.numeric())
  # # print("i ha")

  # extract PPI for Human proteins that are enriched in two pathways
  string.ppi.df = protein.links %>%
    dplyr::filter(combined_score > string_PPI_score_th) %>%
    dplyr::select(-combined_score) %>%
    # dplyr::filter(from %in% allproteins & to %in% allproteins) %>%
    igraph::graph_from_data_frame(directed = F)

  # get network proximity score
network.proximity(net = string.ppi.df, pathway_a_genes = pathway_a,
                  pathway_b_genes = pathway_b, nPermute = nPermute)
}



network.proximity <- function(net,
                              pathway_a_genes,
                              pathway_b_genes,
                              nPermute){
  d <- pathway_a_genes
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])

  t <- pathway_b_genes
  keep <- which(t %in% V(net)$name)
  t <- unique(t[keep])

  d_td <- get_d(shortest.paths(net, v = t, to=d))
  # return(paste0("[Observed distance: ", d_td))

  z <- permuteTest(net, t, d, d_td, nPermute)
  # One-sided (left-tail) hypothesis testing [alternate hypo: observed z-score < mean of random z-scores]
  p <- pnorm(z)

  # return(paste0("[Observed distance: ", d_td,"; z-score: ", z, "; p-value: ", p,"]"))
  return(c(d_td, z, p))
}

get_d <- function(dist_matrix){
  d <- mean(c(apply(dist_matrix,1,min), apply(dist_matrix,2,min)))

  return(d)
}
getRandD <- function(a_degree, net){
  return(base::sample(which(igraph::degree(net) == a_degree),1, replace = F))
}
permuteTest <- function(net, t, d, d_td, N){
  r <- c()
  for (i in 1:N) {
    t_rand <- sapply(as.numeric(igraph::degree(net, v = t)), getRandD,net)
    d_rand <- sapply(as.numeric(igraph::degree(net, v = d)), getRandD,net)
    d_td_rand <- get_d(shortest.paths(net, v = t_rand %>% unique(), to=d_rand %>% unique()))
    r<-c(r,d_td_rand)
  }
  r[!is.finite(r)] <- NA
  m <- mean(r, na.rm = T)
  s <- sd(r, na.rm = T)
  z <- (d_td - m)/s
  return(z)
}

# permuteTest <- function(net, t, d, d_td, N){
#   require(doParallel)
#   require(foreach)
#   cores=detectCores()
#   cl <- makeCluster(cores[1]/2)
#   registerDoParallel(cl)
#
#
#   r <- foreach(i=1:N,
#                .combine='c',
#                .export = c('getRandD','get_d'),
#                .packages = c('igraph','dplyr')) %dopar% {
#                  t_rand <- sapply(as.numeric(igraph::degree(net, v = t)), getRandD, net)
#                  d_rand <- sapply(as.numeric(igraph::degree(net, v = d)), getRandD, net)
#                  d_td_rand <- get_d(igraph::shortest.paths(net,
#                                                            v = t_rand %>% unique(),
#                                                            to = d_rand %>% unique()))
#                  d_td_rand
#                }
#
#   stopCluster(cl)
#
#   r[!is.finite(r)] <- NA
#   m <- mean(r, na.rm = T)
#   s <- sd(r, na.rm = T)
#   # z <- (d_td - m)/(s/sqrt(N))
#   z <- (d_td - m)/s
#
#   return(z)
# }


