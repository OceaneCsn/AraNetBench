#' Draw the result of a network evaluation
#'
#' @param results result returned by the function `evaluate_network` (of type list)
#' @param export weather or not to save the plot
#' @param filename name of the exported plot
#'
#' @return ggplot/ggraph object
#' @export
#'
#' @examples 
draw_evaluated_network <- function(results, export = FALSE, filename = './evaluated_network.png'){
  
  biggest_component_only = TRUE
  graph <- igraph::graph_from_data_frame(results$edges)
  if(biggest_component_only){
    # removes from the plot the small isolated components of the graph
    components <- igraph::clusters(graph, mode="weak")
    biggest_cluster_id <- which.max(components$csize)
    vert_ids <- igraph::V(graph)[components$membership == biggest_cluster_id]
    graph <- igraph::induced_subgraph(graph, vert_ids)
  }
  
  
  layout <- ggraph::create_layout(graph, layout = "stress")
  layout$gene_type <- ifelse(layout$name %in% regulators, "Regulator", "Target")
  layout$gene_type <- ifelse(stringr::str_detect(layout$name, "mean_"), 
                             "Grouped regulators", layout$gene_type)
  
  layout$regulator <- relevel(factor(layout$gene_type), ref = "Target")
  
  if("Grouped regulators" %in% layout$gene_type)
    layout$gene_type <- forcats::fct_relevel(factor(layout$gene_type), 
                                           c("Regulator", "Grouped regulators", "Target"))
  else
    layout$gene_type <- forcats::fct_relevel(factor(layout$gene_type), 
                                             c("Regulator", "Target"))
  
  
  layout$regulator <-
    ifelse(stringr::str_detect(layout$gene_type, "egulator"),
           "Regulator",
           "Target")
  layout$regulator <-
    relevel(factor(layout$regulator), ref = "Target")
  
  colors <- c(
    "Not supported" = "#393939",
    "TF not studied" = "grey",
    "CHIPSeq" = "#0C7C59",
    "DAPSeq" = "#aceca1",
    "TARGET" = "#FAC05E",
    "CHIPSeq+DAPSeq" = "red",
    "DAPSeq+TARGET" = "blue",
    "CHIPSeq+DAPSeq+TARGET" = "green"
  )
  
  val <- ggraph::ggraph(layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        alpha = !stringr::str_detect(type, "TF not studied"),
        col = type
      ),
      width = 0.71,
      show.legend = TRUE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(shape = regulator, size = gene_type,
                   fill = gene_type),
      stroke = 0.235,
      color = "black", show.legend = F,
    ) + ggraph::theme_graph()+
    ggplot2::scale_shape_manual(values = c(21, 22), guide = "none") +
    ggplot2::scale_size_manual(values = c(2, 3, 1.5) * 0.95, guide = "none") +
    ggraph::scale_edge_color_manual(
      name = "Nature of the \nsupporting knwon\ninteraction",
      values = colors
    ) + ggplot2::labs(title = "Network edges colored according to their experimental evidence",
                      subtitle = paste(round(results$tpr*100, 2), "% of the edges (with validation information available) are supported"))+
    ggraph::scale_edge_alpha_manual(values = c(0.51, 1), guide = "none")
  
  if("Grouped regulators" %in% layout$gene_type)
    val <- val + ggplot2::scale_fill_manual(values = c("#68C872", "darkgreen", "grey"), guide = "none")
  else 
    val <- val + ggplot2::scale_fill_manual(values = c("#68C872", "grey"), guide = "none")
  if(export)
    ggpubr::ggexport(val, filename = filename, res = 400, height = 3000, width = 3000)
  return(val)
}


