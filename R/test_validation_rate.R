

#' Test a network against random
#'
#' @param net dataframe with a column 'from' (regulators) and a column 'to'
#' (targets), representing the inferred network of edges to evaluate
#' @param validation type of edge in the validation database that
#' should be considered to defined a true/supported prediction in the
#' evaluation process.
#' The validation type must be a vector of one or more
#' of the following values : CHIPSeq, DAPSeq, Litterature, TARGET
#' @param N Number of shuffled network to estimate the null ditribution
#' of the validation rate. 200 by default.
#'
#' @return A plot containing the empirical p-pvalue 
#' of the permutation-based test
#' @export
#'
#' @examples
#' data("abiotic_stress_Heat_genes_net")
#' set.seed(999)
#' results <- evaluate_network(abiotic_stress_Heat_genes_net)
#' test_validation_rate(abiotic_stress_Heat_genes_net, N = 100)
test_validation_rate <- function(net,
                                 validation = c("CHIPSeq", "DAPSeq", "TARGET"),
                                 N = 200) {
  from_DIANE <- FALSE
  if (any(stringr::str_detect(net$from, 'mean_'))) {
    from_DIANE = TRUE
    grouped_net <- net
    # each interaction links only one gene to one gene
    net <- flatten_edges(net)
  }
  subset <- validated_edges[validated_edges$from %in% net$from &
                              validated_edges$to %in% net$to, ]
  
  observed <- evaluate_network(net, validation, subset)
  tprs <- c()
  for (i in 1:N) {
    shuffled_net <- net
    shuffled_net$to <- sample(shuffled_net$to, replace = FALSE)
    res <-
      evaluate_network(shuffled_net,
                       validation = validation,
                       subset_validated_edges = subset)
    tprs <- c(tprs, res$tpr)
  }
  pval <- sum(tprs > observed$tpr) / N
  zscore <- (observed$tpr - mean(tprs, na.rm = TRUE)) / sd(tprs, na.rm = TRUE)
  distr <- data.frame("Null distribution" = tprs)
  plot <- ggplot2::ggplot(distr, ggplot2::aes(Null.distribution)) +
    ggplot2::geom_density(fill = "#aceca1", alpha = 0.4) +
    ggplot2::geom_vline(
      xintercept = observed$tpr,
      size = 2,
      color = "#0C7C59",
      show.legend = TRUE
    )  +
    ggplot2::geom_label(
      ggplot2::aes(
        x = observed$tpr,
        label = paste("Observed\nvalidation\nrate:", round(observed$tpr, 3)),
        y = N / 25
      ),
      colour = "#0C7C59",
      angle = -90,
      size = 4,
      nudge_x = -0.0055,
      label.size = 0
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(x = Null.distribution, y = 0),
      color = 'black',
      size = 3,
      alpha = 0.2
    ) +
    ggplot2::labs(
      title = paste(
        "Testing the inferred network's validation rate : P value =",
        round(pval, 5),
        " Z-score = ",
        round(zscore, 3)
      ),
      subtitle = paste(
        "Null distribution of validated edges rates computed on",
        N,
        "shuffled networks"
      )
    ) + ggpubr::theme_pubr() +
    ggplot2::xlab("") + ggplot2::ylab("")
  return(plot)
  
}
