##' semantic similarity between two gene clusters
##'
##' given two gene clusters, this function calculates semantic similarity between them.
##'
##' @title clusterSim
##' @param cluster1 a vector of gene IDs
##' @param cluster2 another vector of gene IDs
##' @param organism species of the Entrez gene IDs. If omitted, it is inferred
##' from `ont`.
##' @param ont one of "HDO", "HPO" and "MPO"
##' @param measure One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.
##' @param combine One of "max", "avg", "rcmax", "BMA" methods, for combining
##' @return similarity
##' @importFrom GOSemSim combineScores
##' @export
##' @author Yu Guangchuang
##' @examples
##' \dontrun{
##'	cluster1 <- c("835", "5261","241", "994")
##'	cluster2 <- c("307", "308", "317", "321", "506", "540", "378", "388", "396")
##'	clusterSim(cluster1, cluster2)
##' }
clusterSim <- function(cluster1, 
                       cluster2, 
                       ont = "HDO",
                       organism = NULL,
                       measure="Wang", 
                       combine="BMA") {
    info <- .resolve_ontology_organism(ont, organism)
    ont <- info$ontology
    organism <- info$organism
    .validate_entrez_ids(cluster1, "cluster1")
    .validate_entrez_ids(cluster2, "cluster2")

    do1 <- sapply(cluster1, gene2DO, organism = organism, ont = ont)
    do2 <- sapply(cluster2, gene2DO, organism = organism, ont = ont)

    do1 <- unlist(do1)
    do2 <- unlist(do2)

    res <- doseSim(DOID1 = do1, DOID2 = do2, measure = measure, ont = ont)
    combineScores(res, combine)
}
