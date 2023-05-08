##' measuring similarities between two DO term vectors.
##'
##' provide two term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 DO term, MPO term or HPO term vector
##' @param DOID2 DO term, MPO term or HPO term vector
##' @param ont one of "DO" and "MPO"
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", "Lin", and "TCSS".
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
doseSim <- function(DOID1,
                  DOID2,
                  measure="Wang",
                  ont = "DO") {
    processTCSS <- FALSE
    if (measure == "TCSS") {
        processTCSS <- TRUE
    } 
    ont <- match.arg(ont, c("DO", "MPO", "HPO"))
    if (ont == "DO") {
        scores <- GOSemSim::termSim(DOID1,DOID2, 
            dodata(processTCSS = processTCSS), measure)
    } else if (ont == "MPO") {
        scores <- GOSemSim::termSim(DOID1,DOID2, 
            mpodata(processTCSS = processTCSS), measure)
    } else if (ont == "HPO") {
        scores <- GOSemSim::termSim(DOID1,DOID2, 
            hpodata(processTCSS = processTCSS), measure)        
    }
    if(length(scores) == 1)
        scores <- as.numeric(scores)
    return(scores)
}

##' measuring similarities between two MPO term vectors.
##'
##' provide two DO term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 DO term vector
##' @param DOID2 DO term vector
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", "Lin", and "TCSS".
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
doSim <- function(DOID1,
                   DOID2,
                   measure = "Wang") {
    doseSim(DOID1 = DOID1, DOID2 = DOID2, measure = measure, ont = "DO")                
}

##' measuring similarities between two MPO term vectors.
##'
##' provide two MPO term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 MPO term vector
##' @param DOID2 MPO term vector
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", "Lin", and "TCSS".
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
mpoSim <- function(DOID1,
                   DOID2,
                   measure = "Wang") {
    doseSim(DOID1 = DOID1, DOID2 = DOID2, measure = measure, ont = "MPO")                
}

##' measuring similarities between two MPO term vectors.
##'
##' provide two HPO term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 HPO term vector
##' @param DOID2 HPO term vector
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", "Lin", and "TCSS".
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
hpoSim <- function(DOID1,
                   DOID2,
                   measure = "Wang") {
    doseSim(DOID1 = DOID1, DOID2 = DOID2, measure = measure, ont = "HPO")                
}