##' enrichment map
##'
##' enrichment map
##' @title enrichMap
##' @param x gseaResult or enrichResult object
##' @param fixed if set to FALSE, will invoke tkplot
##' @return figure
##' @importFrom igraph delete.edges
##' @export
##' @author G Yu
enrichMap <- function(x, fixed=TRUE) {
    if (is(x, "gseaResult")) {
        geneSets <- x@geneSets
    }
    if (is(x, "enrichResult")) {
        geneSets <- x@geneInCategory
    }
    y <- summary(x)
    pvalue <- y$pvalue
    
    id <- y[,1]
    geneSets <- geneSets[id]

    n <- nrow(y)
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
    E(g)$width=(wd[,3])*3
    g <- delete.edges(g, E(g)[wd[,3] < 0.2])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
    V(g)$color <- seq_gradient_pal("red", "grey")(pvalue[idx])
    netplot(g, vertex.label.font=1, vertex.label.color="black", fixed=fixed)
}

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}
