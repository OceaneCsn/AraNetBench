#' Deals with grouped nodes of a network
#'
#' The line :
#' from mean_ATGSJHKS-ATGODLKSL to ATGDHJKHK
#'
#' becomes :
#'
#' from ATGSJHKS to ATGDHJKHK
#' from ATGODLKSL to ATGDHJKHK
#'
#' @param net network containing grouped regulators,
#' in the form of a dataframe of edges. It must contain a column
#' named 'from' and a column names 'to'n with the AGI of the genes.
#' The grouped nodes must have the format mean_AGI1-AGI2, such as
#' returned by the inference in DIANE
#'
#' @return dataframe of edges with duplicated edges for
#' grouped regulators
#'
#'
#' @export
flatten_edges <- function(net) {
  # interactions with no grouped nodes
  distinct <- net[!stringr::str_detect(net$from, 'mean_') &
                    !stringr::str_detect(net$to, 'mean_'),
                  c("from", "to")]
  
  # interactions with at least one grouped node
  grouped <- net[stringr::str_detect(net$from, 'mean_') |
                   stringr::str_detect(net$to, 'mean_'),]
  
  if (nrow(grouped) > 0) {
    for (i in 1:nrow(grouped)) {
      tf <- grouped[i, "from"]
      targ <- grouped[i, "to"]
      
      # ungroup the TF(s)
      if (stringr::str_detect(tf, "mean_")) {
        tfs <- unlist(strsplit(stringr::str_remove(tf, "mean_"), '-'))
      }
      else
        tfs <- tf
      
      #ungroup the target(s)
      if (stringr::str_detect(targ, "mean_")) {
        targs <- unlist(strsplit(stringr::str_remove(targ, "mean_"), '-'))
      }
      else
        targs <- targ
      
      # add the new interactions with no grouping
      for (tfi in c(tfs)) {
        for (targi in c(targs)) {
          distinct <- rbind.data.frame(distinct, c(tfi, targi))
        }
      }
    }
  }
  
  return(distinct)
}