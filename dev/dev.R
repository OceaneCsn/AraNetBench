# packages imports
usethis::use_package("stringr")
usethis::use_package("igraph")
usethis::use_package("ggraph")
usethis::use_package("forcats")



usethis::use_vignette("Using-AraNetBench")


devtools::document()
devtools::install()



#devtools::load_all()


