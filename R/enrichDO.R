#' DO Enrichment Analysis 
#'
#' Given a vector of genes, this function will return the enrichment DO
#' categories with FDR control.
#'
#' @rdname enrichDO
#' @inheritParams dose_params
#' @return A \code{enrichResult} instance.
#' @export
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @keywords manip
#' @examples
#'
#'	data(geneList)
#' 	gene = names(geneList)[geneList > 1]
#' 	yy = enrichDO(gene, pvalueCutoff=0.05)
#' 	summary(yy)
#'
enrichDO <- function(gene, ont="HDO",
                     organism = "hsa",
                     pvalueCutoff=0.05,
                     pAdjustMethod="BH",
                     universe,
                     minGSSize = 10,
                     maxGSSize = 500,
                     qvalueCutoff=0.2,
                     readable = FALSE){

    enrichDisease(gene = gene,
                  organism = organism,
                  pvalueCutoff = pvalueCutoff,
                  pAdjustMethod = pAdjustMethod,
                  universe = universe,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  qvalueCutoff = qvalueCutoff,
                  readable = readable,
                  ontology = ont)
}
