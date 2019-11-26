# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' @importFrom magrittr %>%

##-----Organizing scaffold, based on initial scaffold and available data
#' Filter edges from a scaffold based on cells and genes of interest.
#'
#' @param scaffold A dataframe with at least two columns. The two first columns should describe the edges - ie. the first column being the source node and the second the target node.
#' @param cois Cells of interest to filter from the scaffold.
#' @param gois Genes of interest to be filtered.
#' @return A dataframe cointaing the edges in the scaffold that have the cois connected to gois, and edges containig at least one gois.
#'
get_scaffold <- function(scaffold, cois, gois){

  #scaffold is a dataframe with two columns, the first being the source and second the target

  scaffold <- as.data.frame(scaffold)
  scaffold.nodes <- union(scaffold[,1], scaffold[,2])

  ## If genes or cells of interest are specified, scaffold edge must contain at least one item of interest

  ## Filter scaffold and cell and gene list to those of interest
  if ( length(cois)>0 ){
    cells <- dplyr::intersect(scaffold.nodes,cois)
  }else{
    cells <- NA
  }
  if ( length(gois) >0 ){
    genes <- dplyr::intersect(scaffold.nodes,gois)
  }
  if ( length(cois)>0 | length(gois) >0 ){
    scaffold.int.nodes <- c(cells,genes)

    scaffold_cois <- scaffold %>%
      dplyr::filter(From %in% scaffold.int.nodes & To %in% scaffold.int.nodes) #filter edges with both nodes of interest (because cois need to be connected to a gois)

    scaffold_genes <- scaffold %>%
      dplyr::filter(From %in% scaffold.nodes & To %in% scaffold.nodes) %>% #filter edges that we have data on
      dplyr::filter(!(From %in% cells) & !(To %in% cells)) %>% #we don't want edges with cells here
      dplyr::filter(From %in% genes | To %in% genes) #edges with at least one gois

    scaffold <- dplyr::union(scaffold_cois, scaffold_genes)
  }
  return(scaffold)
}



##-----Abundance scores
#' Compute the abundance scores of nodes in a network
#'
#' @param dfn A dataframe with the cell abundance and gene expression values, in a tidy format. Each row of the measurement data should contain, at least,
#' the following columns: a sample ID; node ID; expression value of this node for this sample; group of this sample ID.
#' @param node Name of the column with the nodes names.
#' @param ids Name of the column with the ID of each sample.
#' @param exprv Name of the column with the value of expression data.
#' @param group Name of the column with the group information of the sample.
#' @param cois Cells of interest to filter from the scaffold.
#' @param gois Genes of interest to be filtered.
#' @return A dataframe with the ternary bins for each node and the ratio of samples in each group that map to the two upper tertiles.
#'
compute_abundance <- function(dfn, node, ids, exprv, group, cois, gois){

  #renaming all the columns so they match the naming conventions in the functions
  #rename_cols(dfn, node, "Node")
  #rename_cols(dfn, ids, "SampleID")
  #rename_cols(dfn, exprv, "Value")
  #rename_cols(dfn, group, "Group")
  names(dfn)[match(node, names(dfn))] <- "Node"
  names(dfn)[match(ids, names(dfn))] <- "SampleID"
  names(dfn)[match(exprv, names(dfn))] <- "Value"
  names(dfn)[match(group, names(dfn))] <- "Group"

  tertiles <- getNodeTertiles(dfn)

  df_ternary_full_info <- dfn %>%
    dplyr::group_by(Node) %>% tidyr::nest() %>% ## split by nodes
    dplyr::mutate(Bins=purrr::map2(.x = .data$Node, .y = .data$data, tertiles = tertiles,.f = multiBin)) %>% ## add Bins
    dplyr::mutate(RatioPerGroup=purrr::map(.x = .data$Bins, .f = addPerGroupIncludes))

  df_ternary_full_info

}


##-----Concordance scores
#' Compute the concordance scores of nodes in a network
#'
#' @param scaffold A dataframe with a scaffold interaction in each row. Columns should be named "From" and "To".
#' @param nodes_info A dataframe with the ternary bins for each node and the ratio of samples in each group that map to the two upper tertiles. Use the output dataframe from compute_abundance function.
#' @return A dataframe with the concordance score for each edge in the scaffold.
#'
compute_concordance <- function(scaffold, nodes_info){
  ##filter the scaffold for edges with data for both nodes
  scaffold <- scaffold %>%
    dplyr::filter(From %in% nodes_info$Node & To %in% nodes_info$Node)

  ## Lists ternary for easier referencing
  ternary <- as.list(nodes_info$Bins)

  names(ternary) <- nodes_info$Node ## contains all ternary bins

  scaffold_edge_score_full_info <- scaffold %>%
    dplyr::mutate(PairBin=purrr::map2(.x = From, .y = To, ternary=ternary,  .f=join_binned_pair)) %>% # add table of paired bins per sample
    dplyr::mutate(PairTable=purrr::map(PairBin,tabulate_pair)) %>% ## add table of paired bin frequencies
    dplyr::mutate(ratioScores=purrr::map(PairTable,getGroupRatioScores))# %>% ## add concordance ratio

  scaffold_edge_score <- scaffold_edge_score_full_info %>% dplyr::select(From,To,ratioScores) %>% tidyr::unnest(cols=c(ratioScores))

  predicted_network <- scaffold_edge_score

  predicted_network

}

##-----Consolidate nodes table
#' Consolidate the nodes abundance score in a dataframe ready to be saved to file.
#'
#' @param nodes_info A dataframe with the ternary bins for each node and the ratio of samples in each group that map to the two upper tertiles. Use the output dataframe from compute_abundance function.
#' @return A dataframe with abundance score for all nodes across all groups.
#'
get_abundance_table <- function(nodes_info){

  try(if(!identical(colnames(nodes_info), c("Node", "data", "Bins", "RatioPerGroup")))
    stop ("Names in nodes_info should be 'Node', 'data', 'Bins', 'RatioPerGroup'."))

  nodes_info %>% dplyr::select(Node, RatioPerGroup) %>% tidyr::unnest(c(RatioPerGroup))

}

##-----Filter network based on abundance and concordance scores
#' Filter the edges table to nodes above a given abundance ratio and edges above a given concordance level.
#'
#' @param nodes_scores A dataframe with abundance score for all nodes across all groups. Use the output dataframe from get_abundance_table function.
#' @param edges_scores A dataframe with the concordance score for each edge in the scaffold. Use the output dataframe from compute_concordance function.
#' @param abundance The abundance threshold to filter nodes. Should be between 0 and 1.
#' @param concordance The concordance threshold to filter edges. Only positive values.
#' @return A dataframe with all nodes with UpBinRatio above the abundance threshold and edges with ratioScores above the concordance threshold.
#'

filter_network <- function(nodes_scores, edges_scores, abundance, concordance){

  ab_nodes <- nodes_scores %>%
    dplyr::filter(UpBinRatio > abundance)

  network <- edges_scores %>%
    merge(ab_nodes, by.x = c("From", "Group"), by.y = c("Node", "Group")) %>% #filtering nodes that are abundant
    merge(ab_nodes, by.x = c("To", "Group"), by.y = c("Node", "Group")) %>% #filtering nodes that are abundant
    dplyr::filter(ratioScore > concordance) %>% #filtering concordant edges
    dplyr::select(From, To, Group, ratioScore) %>%
    as.data.frame()

  network
}

