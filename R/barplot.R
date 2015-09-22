##' @title fortify
##' @param model enrichResult object
##' @param data not use here
##' @param showCategory Category numbers to show
##' @param order logical
##' @param drop logical 
##' @param ... additional parameter
##' @importFrom ggplot2 fortify
## @S3method fortify enrichResult
##' @method fortify enrichResult
##' @export
fortify.enrichResult <- function(model, data, showCategory=5, order=FALSE, drop=FALSE, ...) {
    res <- summary(model)
    if (drop) {
        res <- res[res$Count != 0, ]
    }
    if (order) {
        idx <- order(res$Count, decreasing=TRUE)
        res <- res[idx,]
    }
    if ( is.numeric(showCategory) ) {
        if ( showCategory <= nrow(res) ) {
            res <- res[1:showCategory,]
        }
    } else { ## selected categories
        res <- res[res$ID %in% showCategory,]
    }
    res$Description <- factor(res$Description,
                              levels=rev(res$Description))
    
    return(res)
}

##' ggplot theme of DOSE
##'
##' @title theme_dose
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 %+replace%
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @export
##' @param font.size font size
theme_dose <- function(font.size=14) {
    theme_bw() %+replace%
    theme(axis.text.x = element_text(colour = "black",
          size = font.size, vjust =1 ),
          axis.title.x = element_text(colour="black",
          size = font.size),
          axis.text.y = element_text(colour = "black",
          size = font.size, hjust =1 ),
          axis.title.y = element_text(colour="black",
          size = font.size, angle=90)
          )
}


##' @importFrom graphics barplot
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 scale_fill_continuous
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
## @S3method barplot enrichResult
##' @title barplot
##' @param height enrichResult object
##' @param font.size font size
##' @param title plot title
##' @param ... other parameter, ignored
##' @method barplot enrichResult
##' @export
barplot.enrichResult <- function(height, font.size=12, title="", ...) {
    ## use *height* to satisy barplot generic definition
    ## actually here is an enrichResult object.
    x <- height
    Description <- Count <- NULL # to satisfy codetools
    p <- ggplot(x, aes(Description, Count), ... )
    p <- p + geom_bar(stat = "identity") + coord_flip() + theme_dose(font.size)

    if("pvalue" %in% colnames(p$data)) {
        pvalue <- NULL # to satisfy codetools
        p <- p + aes(fill=pvalue) +
            scale_fill_continuous(low="red", high="blue")
    } else {
        p <- p+aes(fill=Description) +
            theme(legend.position="none")
    }
    p <- p + ggtitle(title) + xlab("") + ylab("")
    return(p)
}

