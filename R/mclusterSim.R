##' Pairwise semantic similarity for a list of gene clusters
##'
##'
##' @title mclusterSim
##' @param clusters A list of gene clusters
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @param combine One of "max", "avg", "rcmax", "BMA" methods, for combining semantic similarity scores of multiple DO terms associated with gene/protein.
##' @return similarity matrix
##' @importFrom GOSemSim combineScores
##' @export
##' @author Yu Guangchuang
##' @examples
##'
##'	cluster1 <- c("835", "5261","241")
##'	cluster2 <- c("578","582")
##'	cluster3 <- c("307", "308", "317")
##'	clusters <- list(a=cluster1, b=cluster2, c=cluster3)
##'	mclusterSim(clusters, measure="Wang")
##'
mclusterSim <- function(clusters, measure="Wang", combine="BMA") {
    cluster_dos <- list()
    for (i in seq_along(clusters)) {
        cluster_dos[[i]] <- unlist(sapply(clusters[[i]], gene2DO))
    }
    n <- length(clusters)
    scores <- matrix(NA, nrow=n, ncol=n)
    rownames(scores) <- names(clusters)
    colnames(scores) <- names(clusters)

    for (i in seq_along(cluster_dos)) {
        do1 <- cluster_dos[[i]]
        do1 <- do1[!is.na(do1)]
        for (j in 1:i) {
            do2 <- cluster_dos[[j]]
            do2 <- do2[!is.na(do2)]
            if (length(do1) != 0 && length(do2) != 0) {
                s <- doSim(do1, do2, measure = measure)
                scores[i,j] <- combineScores(s, combine)
                if (i != j) {
                    scores[j, i] <- scores[i, j]
                }
            }
        }
    }
    removeRowNA <- apply(!is.na(scores), 1, sum)>0
    removeColNA <- apply(!is.na(scores), 2, sum)>0
    return(scores[removeRowNA, removeColNA])
}
