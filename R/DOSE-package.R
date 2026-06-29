#' @keywords internal
#' @importClassesFrom enrichit enrichResult gseaResult
#' @importFrom methods is new
#' @importFrom utils stack
"_PACKAGE"

# Declare global variables to suppress R CMD check NOTEs
utils::globalVariables(c(
  "DGN_PATHID2EXTID", "DGN_PATHID2NAME",
  "NCG_PATHID2EXTID", "NCG_PATHID2NAME",
  "VDGN_PATHID2EXTID", "VDGN_PATHID2NAME"
))


#' Datasets
#'
#' Information content and DO term to entrez gene IDs mapping
#'
#'
#' @name DataSet
#' @aliases  geneList NCG_EXTID2PATHID NCG_PATHID2EXTID NCG_PATHID2NAME DGN_EXTID2PATHID DGN_PATHID2EXTID DGN_PATHID2NAME VDGN_EXTID2PATHID VDGN_PATHID2EXTID VDGN_PATHID2NAME
#' @docType data
#' @keywords datasets
NULL
