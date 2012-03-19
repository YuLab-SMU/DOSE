##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 opts
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_text
##' @importFrom ggplot2 theme_bw
.barplot <- function(result,
                    title,
                    font.size=12) {

    Description <- Count <- NULL
    pg <- ggplot(result, aes(x=Description, y=Count)) +
        geom_bar() +
            coord_flip() +
                xlab("") + ylab("") +
                    opts(axis.text.x = theme_text(colour = "black",
                         size = font.size, vjust =1 ),
                         axis.text.y = theme_text(colour = "black",
                         size = font.size, hjust =1 ),
                         title = title) +
                             theme_bw()

    return(pg)
}

##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 scale_fill_continuous
##' @importFrom ggplot2 aes
##' @author Guangchuang Yu \url{http://ygc.name}
barplot.enrichResult <- function(x,
                                 order = FALSE,
                                 drop = FALSE,
                                 showCategory = 5,
                                 title = "",
                                 font.size = 12) {


    res <- summary(x)
    if (drop == TRUE) {
        res <- res[res$Count != 0, ]
    }
    if ( showCategory <= nrow(res) ) {
        res <- res[1:showCategory,]
    }
    if (order == TRUE) {
        idx <- order(res$Count)
        res <- res[idx,]
    }
    res$Description <- factor(res$Description,
                              levels = as.character(res$Description))
    p <- .barplot(res, title, font.size)
    if("pvalue" %in% colnames(res)) {
        pvalue <- NULL # to satisfy codetools
        p <- p + aes(fill=pvalue) +
            scale_fill_continuous(low="red", high="blue")
    }
    return(p)
}

