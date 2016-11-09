##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",
                "For help: https://guangchuangyu.github.io/", pkgname, "\n\n")

  citation <- paste0("If you use ", pkgname, " in published research, please cite:\n",
                     "Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. ",
                     "DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. ",
                     "Bioinformatics 2015, 31(4):608-609", "\n")

  packageStartupMessage(paste0(msg, citation))

	.initial()
}
