library(DIANE)


data("gene_annotations")
data("regulators_per_organism")


regulators <- regulators_per_organism$`Arabidopsis thaliana`

usethis::use_data(regulators, version = 3, overwrite = T)
