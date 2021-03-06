# packages imports
usethis::use_package("stringr")
usethis::use_package("igraph")
usethis::use_package("ggraph")
usethis::use_package("forcats")
usethis::use_package("ggplot2")
usethis::use_package("Cairo")
usethis::use_package("ggpubr")

usethis::use_vignette("Using-AraNetBench")


devtools::document()
devtools::install()

usethis::use_gpl3_license()

pkgdown::build_site()

devtools::build_vignettes()


#devtools::load_all()


