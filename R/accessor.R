##' @method as.data.frame enrichResult
##' @export
as.data.frame.enrichResult <- function(x, ...) {
    as.data.frame(x@result, ...)
}

##' @method as.data.frame gseaResult
##' @export
as.data.frame.gseaResult <- function(x, ...) {
    as.data.frame(x@result, ...)
}

##' @method geneID enrichResult
##' @export
geneID.enrichResult <- function(x) as.character(x@result$geneID)

##' @method geneID gseaResult
##' @export
geneID.gseaResult <- function(x) as.character(x@result$core_enrichment)


##' @method geneInCategory enrichResult
##' @export
##' @importFrom stats setNames
geneInCategory.enrichResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))

##' @method geneInCategory gseaResult
##' @export
geneInCategory.gseaResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))

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
    gc <- geneInCategory(x)
    if (!i %in% names(gc))
        stop("input term not found...")
    gc[[i]]
}


##' @method [[ gseaResult
##' @export
`[[.gseaResult` <- function(x, i) {
    gc <- geneInCategory(x)
    if (!i %in% names(gc))
        stop("input term not found...")
    gc[[i]]
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



