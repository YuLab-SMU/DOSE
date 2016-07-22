##' @method [ enrichResult
##' @export
`[.enrichResult` <- function(x, i, j) {
              x@result[i,j]
}

##' @method [ gseaResult
##' @export
`[.gseaResult` <- `[.enrichResult`


##' @method $ enrichResult
##' @export
`$.enrichResult` <-  function(x, name) {
    x@result[, name]
}

##' @method $ gseaResult
##' @export
`$.gseaResult` <- `$.enrichResult`


##' @method [[ enrichResult
##' @export
`[[.enrichResult` <- function(x, i) {
    if (!i %in% names(x@geneInCategory))
        stop("input term not found...")
    x@geneInCategory[[i]]
}


##' @method [[ gseaResult
##' @export
`[[.gseaResult` <- function(x, i) {
    if (!i %in% names(x@core_enrichment))
        stop("input term not found...")
    x@core_enrichment[[i]]
}


##' @importFrom utils head
##' @method head enrichResult
##' @export
head.enrichResult <- function(x, n=6L, ...) {
    head(x@result, n, ...)
}

##' @method head gseaResult
##' @export
head.gseaResult <- head.enrichResult

##' @importFrom utils tail
##' @method tail enrichResult
##' @export
tail.enrichResult <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}

##' @method tail gseaResult
##' @export
tail.gseaResult <- tail.enrichResult

##' @method dim enrichResult
##' @export
dim.enrichResult <- function(x) {
    dim(x@result)
}

##' @method dim gseaResult
##' @export
dim.gseaResult <- dim.enrichResult



