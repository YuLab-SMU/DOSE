##' upsetplot
##'
##' 
##' @title upsetplot
##' @param x enrichResult object
##' @param n number of categories
##' @return plot 
##' @author Guangchuang Yu
##' @examples
##' \dontrun{
##' require(DOSE)
##' data(geneList)
##' de=names(geneList)[1:100]
##' x <- enrichDO(de)
##' upsetplot(x, 8)
##' }
##' @export
upsetplot <- function(x, n=10) {
    if (! is(x, "enrichResult")) {
        stop("only 'enrichResult' object supported...")
    }
    
    df <- summary(x)
    id <- df$ID[1:n]
    des <- df$Description[1:n]
    glist <- x@geneInCategory[id]
    g <- unique(unlist(glist))


    dat <- matrix(0, nrow=length(g), ncol=length(id))
    rownames(dat) <- g
    for (i in 1:length(id)) {
        dat[glist[[i]], i] <- 1
    }
    colnames(dat) <- des

    ## cols <- ggtree:::color_scale("red", "blue")
    ## pv <- df$pvalue[1:n]
    ## idx <- sapply(pv, function(p) DOSE:::getIdx(p, min(pv), max(pv)))

    ## sets.bar.color = cols[idx],
    UpSetR <- "UpSetR"
    require(UpSetR, character.only = TRUE)
    upset <- eval(parse(text="upset"))
    
    upset(as.data.frame(dat), nsets=n)
}

