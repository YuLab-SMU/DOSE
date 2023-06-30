##' semantic similarity between two gene clusters
##'
##' given two gene clusters, this function calculates semantic similarity between them.
##'
##' @title clusterSim
##' @param cluster1 a vector of gene IDs
##' @param cluster2 another vector of gene IDs
##' @param organism organism
##' @param ont one of "DO" and "MPO"
##' @param measure One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.
##' @param combine One of "max", "avg", "rcmax", "BMA" methods, for combining
##' @return similarity
##' @importFrom GOSemSim combineScores
##' @export
##' @author Yu Guangchuang
##' @examples
##'
##'	cluster1 <- c("835", "5261","241", "994")
##'	cluster2 <- c("307", "308", "317", "321", "506", "540", "378", "388", "396")
##'	clusterSim(cluster1, cluster2)
##'
clusterSim <- function(cluster1, 
                       cluster2, 
                       ont = "DO",
                       organism = "hsa",
                       measure="Wang", 
                       combine="BMA") {
    do1 <- sapply(cluster1, gene2DO, organism = organism)
    do2 <- sapply(cluster2, gene2DO, organism = organism)

    do1 <- unlist(do1)
    do2 <- unlist(do2)

    res <- doseSim(DOID1 = do1, DOID2 = do2, measure = measure, ont = ont)
    combineScores(res, combine)
}
