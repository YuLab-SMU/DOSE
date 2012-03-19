##' measuring similarities between two DO term vectors.
##'
##' provide two DO term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 DO term vector
##' @param DOID2 DO term vector
##' @param method one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @param organism only "human" supported
##' @return score matrix
##' @author Guangchuang Yu \url{http://ygc.name}
doSim <- function(DOID1,
                  DOID2,
                  method="Wang",
                  organism="human") {

    ont <- "DO"
    m <- length(DOID1)
    n <- length(DOID2)
    scores <- matrix(nrow=m, ncol=n)
    rownames(scores) <- DOID1
    colnames(scores) <- DOID2
    for( i in 1:m) {
        for (j in 1:n) {
            if ( is.na(DOID1[i]) || is.na(DOID2[j]) ) {
                scores[i,j] <- NA
            } else {
                if (method == "Wang") {
                    scores[i,j] <- wangMethod(DOID1[i],
                                              DOID2[j],
                                              ont=ont)
                } else {
                    scores[i,j] <- infoContentMethod(DOID1[i],
                                                     DOID2[j],
                                                     ont=ont,
                                                     method=method,
                                                     organism=organism)
                }
            }
        }
    }
    return(scores)
}
