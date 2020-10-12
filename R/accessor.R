##' @method as.data.frame enrichResult
##' @export
as.data.frame.enrichResult <- function(x, ...) {
    x <- get_enriched(x)
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

##' @method geneID compareClusterResult
##' @export
geneID.compareClusterResult <- function(x) as.character(x@compareClusterResult$geneID)


##' @method geneInCategory enrichResult
##' @export
##' @importFrom stats setNames
geneInCategory.enrichResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))

##' @method geneInCategory gseaResult
##' @export
geneInCategory.gseaResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))

##' @method geneInCategory compareClusterResult
##' @export
geneInCategory.compareClusterResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), 
        paste(x@compareClusterResult$Cluster, 
            x@compareClusterResult$ID, sep= "-"))

        
    
    
##' @method [ enrichResult
##' @export
`[.enrichResult` <- function(x, i, j, asis = FALSE, ...) {
    x <- get_enriched(x)
    y <- x@result[i, j, ...]
    if (!asis)
        return(y)
    x@result <- y
    return(x)
}

##' @method [ gseaResult
##' @export
`[.gseaResult` <- function(x, i, j, asis = FALSE, ...) {
    y <- x@result[i, j, ...]
    if (!asis)
        return(y)
    x@result <- y
    return(x)
}


##' @method $ enrichResult
##' @export
`$.enrichResult` <-  function(x, name) {
    x <- get_enriched(x)
    x@result[, name]
}

##' @method $ gseaResult
##' @export
`$.gseaResult` <- function(x, name) {
    x@result[, name]
}



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
    x <- get_enriched(x)
    head(x@result, n, ...)
}

##' @method head gseaResult
##' @export
head.gseaResult <- function(x, n=6L, ...) {
    head(x@result, n, ...)
}

##' @importFrom utils tail
##' @method tail enrichResult
##' @export
tail.enrichResult <- function(x, n=6L, ...) {
    x <- get_enriched(x)
    tail(x@result, n, ...)
}

##' @method tail gseaResult
##' @export
tail.gseaResult <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}

##' @method dim enrichResult
##' @export
dim.enrichResult <- function(x) {
    x <- get_enriched(x)
    dim(x@result)
}

##' @method dim gseaResult
##' @export
dim.gseaResult <- function(x) {
    dim(x@result)
}




##' summary method for \code{gseaResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @return A data frame
##' @exportMethod summary
##' @usage summary(object, ...)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("summary", signature(object="gseaResult"),
          function(object, ...) {
              warning("summary method to convert the object to data.frame is deprecated, please use as.data.frame instead.")
              return(as.data.frame(object, ...))
          }
          )


##' summary method for \code{enrichResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{enrichResult} instance.
##' @param ... additional parameter
##' @return A data frame
##' @exportMethod summary
##' @usage summary(object, ...)
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
setMethod("summary", signature(object="enrichResult"),
          function(object, ...) {
              warning("summary method to convert the object to data.frame is deprecated, please use as.data.frame instead.")
              return(as.data.frame(object, ...))
          }
          )


