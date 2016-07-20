##' convert a list of gene IDs to igraph object.
##'
##'
##' @title convert gene IDs to igraph object
##' @param inputList a list of gene IDs
##' @return a igraph object.
##' @importFrom igraph graph.data.frame
##' @author Guangchuang Yu \url{http://ygc.name}
list2graph <- function(inputList) {
    x <- data.frame()
    for (i in 1:length(inputList)) {
        x=rbind(x,
        data.frame(categoryID=rep(names(inputList[i]),
                   length(inputList[[i]])),
                   Gene=inputList[[i]]))
    }
    g <- graph.data.frame(x, directed=F)
    return(g)
}

##' setting basic attributes of a graph
##'
##' setting size and color of node and edege
##' @title setting.graph.attributes
##' @param g igraph object
##' @param node.size size of node
##' @param node.color color of node
##' @param edege.width edege width
##' @param edege.color color of edege
##' @return igraph object
##' @importFrom igraph V
##' @importFrom igraph "V<-"
##' @importFrom igraph E
##' @importFrom igraph "E<-"
##' @export
##' @author Yu Guangchuang
setting.graph.attributes <- function(g, node.size=8,
                                     node.color="#B3B3B3",
                                     edege.width=2,
                                     edege.color="#8DA0CB") {
    V(g)$size <- node.size
    V(g)$color <- node.color
    V(g)$label <- V(g)$name

    E(g)$width=edege.width
    E(g)$color <- edege.color

    return(g)
}

##' @importFrom scales cscale
##' @importFrom scales seq_gradient_pal
get.col.scale <- function(foldChange, DE.foldChange=FALSE) {
    FC.down <- foldChange[foldChange < 0]
    FC.up <- foldChange[foldChange >=0]

    
    if (DE.foldChange) {
        ## wether foldChange are only DE genes or whole gene sets.
        ## col.down <- cscale(FC.down, seq_gradient_pal("darkgreen", "green"))
        col.down <- cscale(FC.down, seq_gradient_pal("#32FF5C", "#0AFF34"))  
        ## col.up <- cscale(FC.up, seq_gradient_pal("red", "darkred"))
        col.up <- cscale(FC.up, seq_gradient_pal("#FF5C32", "#F52000"))
    } else {
        col.down <- cscale(FC.down, seq_gradient_pal("#32FF5C", "#B3B3B3"))  
        col.up <- cscale(FC.up, seq_gradient_pal("#B3B3B3", "#F52000")) 
    }
    col.scale <- c(col.down, col.up)
    return(col.scale)
}

##' scale color nodes
##'
##' color nodes based on fold change of expression
##' @title scaleNodeColor
##' @param g igraph object
##' @param foldChange fold Change
##' @param node.idx index of node to color
##' @param DE.foldChange logical
##' @importFrom igraph V
##' @importFrom igraph "V<-"
##' @return igraph object
##' @export
##' @author Yu Guangchuang
scaleNodeColor <- function(g, foldChange, node.idx=NULL, DE.foldChange) {
    if (is.null(node.idx)) {
        node.idx <- 1:length(V(g))
    }
    
    gn <- V(g)[node.idx]$name
    if(missing(DE.foldChange) || is.null(DE.foldChange)) {
        if (length(foldChange) > 2*length(gn)) {
            DE.foldChange=FALSE
        } else {
            DE.foldChange=TRUE
        }
    }
    col.scale <- get.col.scale(foldChange, DE.foldChange)
    
    V(g)[node.idx]$color <- col.scale[gn]
    V(g)[node.idx]$color[is.na(V(g)[node.idx]$color)] = "#B3B3B3"
    return(g)
}

##' plot network
##'
##' plot network of igraph object
##' @title netplot
##' @param g igraph object
##' @param vertex.label.font font size
##' @param vertex.label.color font text color
##' @param vertex.label.cex cex of vertex label
##' @param layout layout
##' @param foldChange fold change
##' @param fixed logical
##' @param col.bin number of legend color bin
##' @param legend.x x-axis position of legend
##' @param legend.y y-axis position of legend
##' @importFrom igraph tkplot
##' @importFrom igraph plot.igraph
##' @importFrom igraph layout.fruchterman.reingold
##' @importFrom graphics hist
##' @importFrom graphics points
##' @importFrom graphics text
##' @export
##' @return plot
##' @author Yu Guangchuang
netplot <- function(g,
                    vertex.label.font=2,
                    vertex.label.color='#666666',
                    vertex.label.cex=1.5,
                    layout=layout.fruchterman.reingold,
                    foldChange=NULL,
                    fixed=TRUE,
                    col.bin=10,
                    legend.x=1,
                    legend.y=1) {
    if (fixed){
        plot.igraph(g,
                    vertex.label.font=vertex.label.font,
                    vertex.label.color=vertex.label.color,
                    vertex.label.cex=vertex.label.cex,
                    vertex.frame.color=V(g)$color,
                    layout=layout)
        ## add legend
        if (!is.null(foldChange)) {
            ## gn <- V(g)$name
            ## fc <- foldChange[gn]
            ## fc <- fc[!is.na(fc)]
            fc <- foldChange
            lbs <- hist(fc, breaks=col.bin-1, plot=FALSE)$breaks
            col.legend <- get.col.scale(lbs)

            x <- seq(from=legend.x, by=0.03, length.out=length(col.legend))
            y <- rep(legend.y, length(col.legend))
            points(x, y, pch=15, col=col.legend, cex=2)

            idx <- c(1, seq(4, length(col.legend)-1, by=3), length(col.legend))
            text(x=x[idx],
                 y=rep(legend.y-0.05, length(idx)),
                 label=lbs[idx],
                 cex = 0.8)

            text(x=mean(x),
                 y=legend.y+0.05,
                 labels="Fold Change",
                 cex=0.8, font=2)
        }
    } else {
        tkplot(g,
               vertex.label.font=vertex.label.font,
               vertex.label.color=vertex.label.color,
               vertex.label.cex=vertex.label.cex,
               vertex.frame.color=V(g)$color,
               layout=layout)
    }
    invisible(g)
}

##' plot function of gene Concept Net.
##'
##'
##' @title cnetplot_internal
##' @param inputList a list of gene IDs
##' @param categorySize setting category size
##' @param showCategory number of categories to plot
##' @param pvalue pvalue
##' @param foldChange  fold Change
##' @param fixed logical
##' @param DE.foldChange logical
##' @param ... additional parameters
##' @return plotted igraph object.
##' @importFrom igraph V
##' @importFrom igraph "V<-"
##' @importFrom igraph degree
##' @author Guangchuang Yu \url{http://ygc.name}
cnetplot_internal <- function(inputList,
                              categorySize="geneNum",
                              showCategory=5,
                              pvalue=NULL,
                              foldChange=NULL,
                              fixed=TRUE, 
                              DE.foldChange = NULL, ...) {

    if (is.numeric(showCategory)) {
        inputList <- inputList[1:showCategory]
        if (!is.null(pvalue)) {
            pvalue <- pvalue[1:showCategory]
        }
    } else {## selected categories
        inputList <- inputList[showCategory]
        if ( !is.null(pvalue) ) {
            pvalue <- pvalue[showCategory]
        }
    }

    ## generate graph object
    g <- list2graph(inputList)

    ## setting some attributes
    g <- setting.graph.attributes(g)

    lengthOfCategory <- length(inputList)

    ## scale node colors based on fold change
    if ( !is.null(foldChange) ) {
        node.idx <- (lengthOfCategory+1):length(V(g))
        g <- scaleNodeColor(g, foldChange, node.idx, DE.foldChange)
    }

    ## attributes of Category Node
    V(g)[1:lengthOfCategory]$size=30  ## setting by default.
    V(g)[1:lengthOfCategory]$color= "#E5C494"

    if(is.numeric(categorySize)) {
        V(g)[1:lengthOfCategory]$size=categorySize
    } else {
        if (categorySize == "geneNum") {
            n <- degree(g)[1:lengthOfCategory]
            V(g)[1:lengthOfCategory]$size <- n/sum(n) * 100
        }
        if (categorySize == "pvalue") {
            if (is.null(pvalue)) {
                stop("pvalue must not be null...")
            }
            pScore<- -log10(pvalue)
            V(g)[1:lengthOfCategory]$size <- pScore/sum(pScore) * 100
        }
    }

    netplot(g=g,foldChange=foldChange, fixed=fixed, ...)
}


cnetplot.enrichResult <- function(x,
                                  showCategory = 5,
                                  categorySize = "geneNum",
                                  foldChange   = NULL,
                                  fixed        = TRUE,
                                  ...) {
    res <- summary(x)
    if (is(x, "enrichResult")) {
        gc <- x@geneInCategory
    } else if (is(x, "gseaResult")) {
        gc <- x@core_enrichment
    } else {
        stop("x should be an 'enrichResult' or 'gseaResult' object...")
    }
    
    if ("pvalue" %in% names(res)) {
        y <- res[res$ID %in% names(gc),
                 c("ID", "Description", "pvalue")]
    } else {
        y <- res[res$ID %in% names(gc),
                 c("ID", "Description")]
        
    }
    
    gc <- gc[as.character(y$ID)]
    names(gc) <- y$Description
    
    if ( is.numeric(showCategory) && (showCategory > length(gc)) ) {
        showCategory = length(gc)
    }
    
    if (categorySize == "pvalue") {
        pvalue <- y$pvalue
        names(pvalue) <- y$Description
    } else {
        pvalue <- NULL
    }

    readable <- x@readable
    organism <- x@organism
    if (readable & (!is.null(foldChange) ) ){
        gid <- names(foldChange)
        if (is(x, 'gseaResult')) {
            ii <- gid %in% names(x@geneList)
        } else {
            ii <- gid %in% x@gene
        }
        gid[ii] <- x@gene2Symbol[gid[ii]]
        names(foldChange) <- gid
    }
    
    cnetplot_internal(inputList=gc,
                      showCategory=showCategory,
                      categorySize=categorySize,
                      pvalue=pvalue,
                      foldChange=foldChange,
                      fixed=fixed, ...)
}
