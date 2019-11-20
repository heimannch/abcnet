#------utils functions

## this assert function should maybe go in functions/utils.R
assert_list_has_names <- function(lst, names){
  missing_names <- names[!names %in% names(lst)]
  if(length(missing_names) != 0){
    stop("Names/labels not present in the list object: ",
         stringr::str_c(missing_names, collapse = ", "))
  }
}

rename_cols <- function(df,x,y){
  names(df)[match(x, names(df))] <- y
  return(df)
}

##-----Functions for Indiviual Nodes

## helper functions for univariate tertile frequency (should one add asserts to these?)
unibin <- function(df){table(df$Bin)} ## tabulate bin frequencies in Bin column of df
upbinfrac  <- function(v){(as.numeric(v[2]+v[3])/sum(v))} ## for 3-vector, ratio of upper two values to all
#keep.node <- function(val,upbinfrac.threshold){val>=upbinfrac.threshold} ## helper function for thresholding upper bin ratio

## Get tertiles for each node
getNodeTertiles <- function(df){
  quantz <- df %>% group_by(Node) %>% summarise(q0=quantile(Value,0,na.rm=T),
                                                q33=quantile(Value,0.33,na.rm=T),
                                                q66=quantile(Value,0.66,na.rm=T),
                                                q1=quantile(Value,1,na.rm=T))
  tertiles <- as.matrix(quantz[,2:5]) ; rownames(tertiles) <- quantz$Node
  tertiles
}

## Bin values according to tertiles
##
getbin <- function(value,node_label, tertiles){
  cut(value,tertiles[node_label,])
}

## Bin all node values for a particular node
## input
## s: node name
## df: data frame of ParticipantBarcode,Group,Value
## requires tertiles for getbin function
multiBin <- function(s,df,tertiles, byTwoGroups = FALSE){

  vals <- df$Value
  beans <- getbin(vals,s, tertiles)

  tibble(SampleID=df$SampleID,
         Group=df$Group,Bin=beans)

}

## Calculate node inclusion, per Group
## Acts on a df ParticipantBarcode,Group,Bin
## Intermediate steps calculate bin fractions
## Outputs df with Group, and Include, where latter is logical
addPerGroupIncludes <- function(df){

    ratios <- df %>% group_by(Group) %>% tidyr::nest() %>% ## split by Group, now have tibble with one row per Group; tibble named data; has ParticipantBarcode, Bin
      mutate(BinDistribution=purrr::map(data,unibin)) %>% ## add BinDistribution tibble with frequency of samples in each tertile bin
      mutate(UpBinRatio=purrr::map_dbl(BinDistribution,upbinfrac)) %>% ## add UpBinRatio with ratio of upper two bins to all (three) bins
      select(Group,UpBinRatio) %>% arrange(Group)

  return(ratios)
}

## Functions for Node Pairs

## for two nodes create data frame of bins for both pair members
## requires object ternary
join_binned_pair <- function(
  NodeA,
  NodeB,
  ternary,
  group_column = "Group",
  id_column = "SampleID"){

  assert_list_has_names(ternary,c(NodeA,NodeB))

  dplyr::inner_join(ternary[[NodeA]],ternary[[NodeB]],by=c(id_column,group_column)) %>%
    dplyr::select(-id_column)
}

tabulate_pair <- function (df){

    df %>% group_by(Group) %>% tidyr::nest() %>% mutate(t=purrr::map(data,table)) %>% select(Group,t)

}

## diagonal over off-diagonal cross-tabulated bin counts
## uses pseudocount for denominator
getRatioScore <- function(Xtab,pseudocount=1){
  ratioScore <- (Xtab[3,3] + Xtab[1,1]) / (Xtab[3,1] + Xtab[1,3]+pseudocount)
  # if (ratioScore < 1) {
  #   ratioScore <- -(1/ratioScore)
  # }
  ratioScore
}

getGroupRatioScores <- function(pt){

    pt %>% ungroup() %>% mutate(ratioScore=purrr::map(t,getRatioScore)) %>%
      transmute(Group,ratioScore=as.numeric(ratioScore)) %>%
      arrange(Group)

}





#---------------------------
# Functions used in iAtlas
#---------------------------

## get nodes of a particular type
## requires scaffold
# get_nodes_of_type <- function(type_of_node){
#   unique(sort(c(names(which(node_type_of_node[scaffold$From]==type_of_node)),
#                 names(which(node_type_of_node[scaffold$To]==type_of_node)))))
# }
##-----Organizing cell and gene expression data

# get_participants <- function(group_df,
#                              group_col,
#                              byImmune = FALSE){
#   # The participant list is obtained from the data frame that has the groups
#
#   participants <- group_df$ParticipantBarcode
#   groups <- group_df[[group_col]]
#   group_of_participant <- groups ; names(group_of_participant) <- participants
#
#   return(group_of_participant)
# }

#organizing the TCGA data in the right format
# get_dfn <- function(dfc_in,
#                     dfe_in,
#                     group_of_participant,
#                     cells,
#                     genes,
#                     group_df,
#                     byImmune = FALSE){
#
#   set.seed(42)
#
#   participants <- names(group_of_participant)
#
#   dfc <- dfc_in %>% filter(ParticipantBarcode %in% participants) %>%
#     select(ParticipantBarcode,paste(cells,".Aggregate2",sep=""))
#   colnames(dfc) <- gsub(".Aggregate2","",colnames(dfc))
#   dfc <- dfc %>% dplyr::mutate(Group=group_of_participant[ParticipantBarcode])
#
#   if(byImmune == FALSE){
#
#     dfclong <- dfc %>% tidyr::gather(Cell,Cell_Fraction,-c(ParticipantBarcode,Group))
#
#   }else{
#     immune_group <- get_participants(group_df, "Subtype_Immune_Model_Based")
#     dfc <- dfc %>% dplyr::mutate(Immune=immune_group[ParticipantBarcode])
#
#     dfclong <- dfc %>% tidyr::gather(Cell,Cell_Fraction,-c(ParticipantBarcode,Group, Immune))
#   }
#
#   dfclong <- dfclong %>% dplyr::mutate(Cell_Fraction=Cell_Fraction+rnorm(mean=0, sd=0.0001,nrow(.)))
#
#   dfclong.generic <- dfclong %>% rename(Node=Cell,Value=Cell_Fraction)
#
#   #computing gene expression data
#
#   dfe <- dfe_in %>% filter(ParticipantBarcode %in% participants) %>% select(ParticipantBarcode,as.vector(genes))
#   dfelong <- dfe %>% tidyr::gather(Gene,Expression,-ParticipantBarcode)
#   dfelong <- dfelong %>% dplyr::mutate(ExpLog2 = log2(Expression+1)+rnorm(mean=0, sd=0.0001,nrow(.))) %>%
#     dplyr::select(ParticipantBarcode,Gene,Expression=ExpLog2) %>%
#     dplyr::mutate(Group=group_of_participant[ParticipantBarcode])
#
#   if(byImmune == FALSE){
#     dfelong <- dfelong %>%
#       dplyr::select(ParticipantBarcode,Group,Gene,Expression)
#   }else{
#     immune_group <- get_participants(group_df, "Subtype_Immune_Model_Based")
#
#     dfelong <- dfelong %>%
#       dplyr::mutate(Immune = immune_group[ParticipantBarcode]) %>%
#       dplyr::select(ParticipantBarcode,Group, Immune,Gene,Expression)
#   }
#
#   dfelong.generic <- dfelong %>% rename(Node=Gene,Value=Expression)
#   dfelong.generic
#
#   dfn <- dplyr::bind_rows(dfelong.generic, dfclong.generic)
#
#   return(dfn)
# }
#
## Retain concordance ratio if ratio meets concordance threshold and both pair members meet include criterion
#
# filterRatioWithIncludes <- function (rscores,from_node,to_node,feature_include,concordance_treshold){ # concordance_treshold=1.62
#   from_include <- feature_include[[from_node]]
#   to_include <- feature_include[[to_node]]
#   b <- inner_join(from_include,to_include,by="Group") %>% mutate(Both=Include.x & Include.y) %>%
#     select(Group,Both)
#   g <- function(val,logval){if(logval){val}else{NA}} ## helper for combining two includes
#   h <- function(val){if (is.na(val)){NA}else if(val>concordance_treshold){val}else{NA}} # helper for threshold filter
#   inner_join(rscores,b,by="Group") %>%
#     mutate(RatioFiltered=purrr::map2_dbl(ratioScore,Both,g)) %>%
#     mutate(RatioFiltered=purrr::map_dbl(RatioFiltered,h)) %>%
#     select(Group,RatioFiltered)
# }


