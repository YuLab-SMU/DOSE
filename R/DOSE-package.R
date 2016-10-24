#' Disease Ontology Semantic and Enrichment analysis
#' Implemented five methods proposed by Resnik, Schlicker, Jiang, Lin and Wang respectively for measuring DO semantic similarities, and hypergeometric test for enrichment analysis.
#'
#' This package is designed to estimate DO-based semantic similarity
#' measurement and enrichment analysis.
#'
#' \tabular{ll}{ Package: \tab DOSE\cr Type: \tab Package\cr Version: \tab
#' 2.3.5\cr Date: \tab 2-27-2012\cr biocViews:\tab Bioinformatics,
#' Annotation\cr Depends:\tab \cr Imports: \tab methods, AnnotationDbi,
#' DO.db\cr Suggests:\tab clusterProfiler, GOSemSim\cr License: \tab
#' Artistic-2.0\cr }
#'
#' @name DOSE-package
#' @aliases DOSE-package DOSE
#' @docType package
#' @author Guangchuang Yu, Li-Gen Wang
#'
#' Maintainer: Guangchuang Yu \email{guangchuangyu@@gmail.com}
#' @seealso \linkS4class{enrichResult}
#' @keywords package
NULL


#' Datasets
#'
#' Information content and DO term to entrez gene IDs mapping
#'
#'
#' @name DataSet
#' @aliases EG2DO DO2EG EG2ALLDO DO2ALLEG DOIC dotbl
#' EG2DOLite DOLite2EG DOLiteTerm geneList
#' NCG_EXTID2PATHID NCG_PATHID2EXTID NCG_PATHID2NAME
#' DGN_EXTID2PATHID DGN_PATHID2EXTID DGN_PATHID2NAME
#' VDGN_EXTID2PATHID VDGN_PATHID2EXTID VDGN_PATHID2NAME
#' @docType data
#' @keywords datasets
NULL
