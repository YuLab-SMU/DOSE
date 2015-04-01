##' enrichment map
##'
##' enrichment map
##' @title enrichMap
##' @param x gseaResult or enrichResult object
##' @param n maximum number of category to shown
##' @param fixed if set to FALSE, will invoke tkplot
##' @param ... additional parameter
##' @return figure
##' @importFrom igraph delete.edges
##' @export
##' @author G Yu
enrichMap <- function(x, n = 50, fixed=TRUE, ...) {
    if (is(x, "gseaResult")) {
        geneSets <- x@geneSets
    }
    if (is(x, "enrichResult")) {
        geneSets <- x@geneInCategory
    }
    y <- summary(x)
    if (nrow(y) < n) {
        n <- nrow(y)
    }
    y <- y[1:n,]
    
    pvalue <- y$pvalue
    
    id <- y[,1]
    geneSets <- geneSets[id]
 
    n <- nrow(y) #
    w <- matrix(NA, nrow=n, ncol=n)
    colnames(w) <- rownames(w) <- y$Description
    
    for (i in 1:n) {
        for (j in i:n) {
            w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
        }
    }
   
    wd <- melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    wd <- wd[!is.na(wd[,3]),]
    g <- graph.data.frame(wd[,-3], directed=F)
    E(g)$width=sqrt(wd[,3]*20)
    g <- delete.edges(g, E(g)[wd[,3] < 0.2])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))

    cols <- color_scale("red", "#E5C494")
    
    V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))]
    ## seq_gradient_pal("red", "grey")(pvalue[idx])
    netplot(g, vertex.label.font=1, vertex.label.color="black", fixed=fixed, ...)
}

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

color_scale <- function(c1="grey", c2="red") {
    pal <- colorRampPalette(c(c1, c2))
    colors <- pal(100)
    return(colors)
}

getIdx <- function(v, MIN, MAX) {
    if ( MIN == MAX ) {
        return(100)
    }
    intervals <- seq(MIN, MAX, length.out=100)
    max(which(intervals <= v))
}
