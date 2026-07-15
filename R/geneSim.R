##' measuring similarities bewteen two gene vectors.
##'
##' provide two entrez gene vectors, this function will calculate their similarity.
##' @title geneSim
##' @param geneID1 entrez gene vector
##' @param geneID2 entrez gene vector
##' @param organism species of the Entrez gene IDs. If omitted, it is inferred
##' from `ont`.
##' @param ont one of "HDO" and "MPO"
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @param combine One of "max", "avg", "rcmax", "BMA" methods, for combining semantic similarity scores of multiple DO terms associated with gene/protein.
##' @return score matrix
##' @importFrom GOSemSim combineScores
##' @export
##' @examples
##' g <- c("835", "5261","241", "994")
##' geneSim(g)
##' @author Guangchuang Yu \url{https://yulab-smu.top}
geneSim <- function(geneID1,
                    geneID2=NULL,
                    ont = "HDO",
                    organism = NULL,
                    measure="Wang",
                    combine="BMA") {

    info <- .resolve_ontology_organism(ont, organism)
    ont <- info$ontology
    organism <- info$organism
    .validate_entrez_ids(geneID1, "geneID1")
    if (!is.null(geneID2)) .validate_entrez_ids(geneID2, "geneID2")

    DOID1 <- lapply(geneID1, gene2DO, organism = organism, ont = ont)
    if (is.null(geneID2)) {
        geneID2 <- geneID1
        DOID2 <- DOID1
    } else {
        DOID2 <- lapply(geneID2, gene2DO, organism = organism, ont = ont)
    }

    m <- length(geneID1)
    n <- length(geneID2)
    scores <- matrix(NA, nrow=m, ncol=n)
    rownames(scores) <- geneID1
    colnames(scores) <- geneID2

    for (i in 1:m) {
        if (length(geneID1) == length(geneID2) && all(geneID1 == geneID2)) {
           nn <- i
           flag <- TRUE
        } else {
            flag <- FALSE
            nn <- n
        }

        for (j in 1:nn) {
            if(any(!is.na(DOID1[[i]])) &&  any(!is.na(DOID2[[j]]))) {
                s <- doseSim(DOID1[[i]],
                           DOID2[[j]],
                           measure = measure,
                           ont = ont
                           )
                scores[i,j] = combineScores(s, combine)
                if (flag == TRUE && j != i) {
                    scores[j, i] <- scores[i,j]
                }
            }
        }
    }
    if (nrow(scores) == 1 & ncol(scores) == 1)
        scores = as.numeric(scores)
    return(scores)
}
