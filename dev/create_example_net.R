# network from current ongoing work on CO2*N
# load("D:/These/Combi/Results/network_N_CO2.RData")
# 
# network_N_CO2$network_data$edges

library(DIANE)
data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")

# getting raw counts from DIANE's demo data
tcc_object <- list(counts = abiotic_stresses$raw_counts)
threshold = 10*length(abiotic_stresses$conditions)
tcc_object$counts <- tcc_object$counts[rowSums(tcc_object$counts) > threshold,]
normalized_counts <- tcc_object$counts

fit <- DIANE::estimateDispersion(tcc = tcc_object)


topTags <- DIANE::estimateDEGs(fit, reference = "M", 
                               perturbation = "MH", 
                               p.value = 0.05, lfc = 1.5)

genes <- get_locus(topTags$table$genes)
regressors <- intersect(genes, regulators_per_organism[[
  "Arabidopsis thaliana"]])

aggregated_data <- aggregate_splice_variants(data = normalized_counts)

grouping <- DIANE::group_regressors(aggregated_data, genes, 
                                    regressors, corr_thr = 0.9)

grouped_counts <- grouping$counts
grouped_targets <- grouping$grouped_genes
grouped_regressors <- grouping$grouped_regressors


mat <- DIANE::network_inference(grouped_counts, 
                                conds = abiotic_stresses$conditions, 
                                targets = grouped_targets, 
                                regressors = grouped_regressors, 
                                importance_metric = "MSEincrease_oob",
                                nCores = 5, verbose = FALSE, nTrees = 2000)


network <- DIANE::network_thresholding(mat, n_edges = length(genes))
data <- network_data(network, regulators_per_organism[["Arabidopsis thaliana"]], gene_annotations$`Arabidopsis thaliana`)
net <- data$edges[,c("from", "to")]
abiotic_stress_Heat_genes_net <- net

usethis::use_data(abiotic_stress_Heat_genes_net, version = 3, overwrite = T)


# save(net, file = "abiotic_stress_Heat_genes_network.rdata")
a <- draw_network(data$nodes, data$edges);a
