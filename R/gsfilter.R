##' filter enriched result by gene set size or gene count
##'
##' 
##' @title gsfilter
##' @param x instance of enrichResult or compareClusterResult
##' @param by one of 'GSSize' or 'Count'
##' @param min minimal size
##' @param max maximal size
##' @return update object
##' @export
##' @author Guangchuang Yu
gsfilter <- function(x, by="GSSize", min=NA, max=NA) {
    by <- match.arg(by, c("GSSize", "Count"))
    if (is(x, "enrichResult")) {
        result <- x@result
    } else {
        result <- x@compareClusterResult
    }
    
    if (by == "GSSize") {
        n <- as.numeric(gsub("^(\\d+)/\\d+$",'\\1',  as.character(result$BgRatio)))
    } else {
        n <- result$Count
    }
    if (!is.na(min)) {
        min_lidx <- n >= min
    } else {
        min_lidx <- TRUE
    }
    
    if (!is.na(max)) {
        max_lidx <- n <= max
    } else {
        max_lidx <- TRUE
    }

    idx <- min_lidx & max_lidx

    if (is(x, "enrichResult")) {
        x@result <- result[idx,]
    } else {
        x@compareClusterResult <- result[idx,]
    }
    return(x)
}
