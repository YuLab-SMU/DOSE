##' semantic similarity between two gene clusters
##'
##' given two gene clusters, this function calculates semantic similarity between them.
##' 
##' @title clusterSim 
##' @param cluster1 a vector of gene IDs
##' @param cluster2 another vector of gene IDs
##' @param measure One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.
##' @param combine One of "max", "average", "rcmax", "BMA" methods, for combining
##' @return similarity
##' @importFrom GOSemSim combineScores
##' @export
##' @author Yu Guangchuang
##' @examples
##'
##'	## cluster1 <- c("835", "5261","241", "994")
##'	## cluster2 <- c("307", "308", "317", "321", "506", "540", "378", "388", "396")
##'	## clusterSim(cluster1, cluster2, ont="MF", organism="human", measure="Wang")
##'
clusterSim <- function(cluster1, cluster2, measure="Wang", combine="BMA") {
    do1 <- sapply(cluster1, gene2DO)
    do2 <- sapply(cluster2, gene2DO)

    do1 <- unlist(do1)
    do2 <- unlist(do2)

    res <- doSim(do1, do2, measure = measure)
    combineScores(res, combine)
}
