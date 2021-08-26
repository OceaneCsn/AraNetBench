#' Evaluate Network
#' 
#' Evaluates an inferred network against validated regulatory interactions
#'
#' @param net dataframe with a column 'from' (regulators) and a column 'to'
#' (targets), representing the inferred network of edges to evaluate
#' @param validation type of edge in the validation database that 
#' should be considered to defined a true/supported prediction in the 
#' evaluation process.
#' The validation type must be a vector of one or more 
#' of the following values : CHIPSeq, DAPSeq, Litterature, TARGET
#' @param subset_validated_edges 
#'
#' @return a list containing true positives, true positive rate, 
#' false positives, false positive rate
#' @export
evaluate_network <-
  function(net,
           validation = c("CHIPSeq", "DAPSeq", "Litterature", "TARGET"), 
           subset_validated_edges = NULL) {
    
    if(!is.null(subset_validated_edges))
      validated_edges <- subset_validated_edges
    
    # ----------------------------- CHECKS ---------------------------------- #
    
    if(any(!validation %in% c("CHIPSeq", "DAPSeq", "Litterature", "TARGET")))
      stop("The validation type must be a vector of one or more of the following values :
           CHIPSeq, DAPSeq, Litterature, TARGET")
    
    if (!("from" %in% colnames(net) &
          "to" %in% colnames(net)))
      stop("The network dataframe should have two columns named 'from' and 'to")
    
    
    # regex check for Arabidopsis AGIs
    matched <-
      sum(stringr::str_detect(
        c(net$from, net$to), pattern = "^AT[[:alnum:]]G[[:digit:]]{5}"))
    
    if (matched != 2 * nrow(net)) {
      if (matched > 0)
        stop("Some of the gene IDs do not match the expected regex for Arabidopsis AGIs")
      else
        stop("None of the gene IDs match the expected regex Arabidopsis AGIs")
    }
    
    # ------------------------------------------------------------------------ #
    
    validated_edges_specific <- 
      validated_edges[validated_edges$type %in% validation, ]
    
    val <-
      merge(net, validated_edges_specific,
            by = c('from', 'to'))
    
    
    # TFs for which we possess validation data
    studied_tfs <- unique(validated_edges_specific$from)
    n_studied_interactions <- sum(net$from %in% studied_tfs)
    
    if(length(studied_tfs) == 0){
      warning("No regulator present in your network was studied in the database \n")
      return(list(
        tp = NA,
        fp = NA,
        tpr = NA,
        fpr = NA
      ))
    }
    
    if(nrow(val) == 0){
      warning("No predicted edge was found in the validation databse... 
              Some coffee to cheer you up? \n")
      return(list(
        tp = 0,
        fp = n_studied_interactions,
        tpr = 0,
        fpr = 1
      ))
    }
    
    val_unique <-
      aggregate(. ~ from + to,
                data = val,
                FUN = paste0,
                collapse = '+')
  
    # true positives
    tp <- nrow(val_unique)
    
    # true positive rate
    if(nrow(val) == 0){
      tpr = 0
    }
    else{
      tpr <- tp / n_studied_interactions
    }
      
    # false positives
    fp <- n_studied_interactions - tp
    
    # false positive rate
    fpr <- fp / n_studied_interactions
    
    results <- list(
      tp = tp,
      fp = fp,
      tpr = tpr,
      fpr = fpr
    )
    return(results)
  }
