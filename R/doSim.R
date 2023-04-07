##' measuring similarities between two DO term vectors.
##'
##' provide two DO term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 DO term or MPO term vector
##' @param DOID2 DO term or MPO term vector
##' @param ont one of "DO" and "MPO"
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", "Lin", and "TCSS".
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
doSim <- function(DOID1,
                  DOID2,
                  measure="Wang",
                  ont = "DO") {
    processTCSS <- FALSE
    if (measure == "TCSS") {
        processTCSS <- TRUE
    } 
    if (ont == "DO") {
        scores <- GOSemSim::termSim(DOID1,DOID2, 
            dodata(processTCSS = processTCSS), measure)
    } else {
        scores <- GOSemSim::termSim(DOID1,DOID2, 
            mpodata(processTCSS = processTCSS), measure)
    }
    if(length(scores) == 1)
        scores <- as.numeric(scores)
    return(scores)
}
