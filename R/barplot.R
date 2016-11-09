##' convert enrichResult object for ggplot2
##'
##'
##' @title fortify
##' @param model enrichResult object
##' @param data not use here
##' @param showCategory Category numbers to show
##' @param by one of Count and GeneRatio
##' @param order logical
##' @param drop logical
##' @param category separate result by 'category' variable
##' @param ... additional parameter
##' @importFrom ggplot2 fortify
## @S3method fortify enrichResult
##' @method fortify enrichResult
##' @export
fortify.enrichResult <- function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, category=NULL, ...) {
    res <- as.data.frame(model)
    if (drop) {
        res <- res[res$Count != 0, ]
    }
    res$GeneRatio <- parse_ratio(res$GeneRatio)

    if (order) {
        if (by == "Count") {
            idx <- order(res$Count, decreasing=TRUE)
        } else {
            idx <- order(res$GeneRatio, decreasing=TRUE)
        }
        res <- res[idx,]
    }

    topN <- function(res, showCategory) {
        if ( is.numeric(showCategory) ) {
            if ( showCategory <= nrow(res) ) {
                res <- res[1:showCategory,]
            }
        } else { ## selected categories
            res <- res[res$ID %in% showCategory,]
        }
        return(res)
    }

    if (is.null(category)) {
        res <- topN(res, showCategory)
    } else {
        lres <- split(res, as.character(res[, category]))
        lres <- lapply(lres, topN, showCategory = showCategory)
        res <- do.call('rbind', lres)
    }

    res$Description <- factor(res$Description,
                              levels=rev(res$Description))

    return(res)
}

##' ggplot theme of DOSE
##'
##' @title theme_dose
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @export
##' @param font.size font size
theme_dose <- function(font.size=14) {
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black",
                                     size = font.size, vjust =1 ),
          axis.text.y = element_text(colour = "black",
                                     size = font.size, hjust =1 ),
          axis.title = element_text(margin=margin(10, 5, 0, 0),
                                    color = "black",
                                    size = font.size),
          axis.title.y = element_text(angle=90)
          )
}


##' barplot of enrichResult
##'
##'
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
##' @param x one of 'Count' and 'GeneRatio'
##' @param colorBy one of 'pvalue', 'p.adjust', 'qvalue'
##' @param showCategory number of categories to show
##' @param font.size font size
##' @param title plot title
##' @param ... other parameter, ignored
##' @method barplot enrichResult
##' @export
barplot.enrichResult <- function(height, x="Count", colorBy='pvalue', showCategory=5, font.size=12, title="", ...) {
    ## use *height* to satisy barplot generic definition
    ## actually here is an enrichResult object.
    object <- height

    colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue"))
    if (x == "geneRatio" || x == "GeneRatio") {
        x <- "GeneRatio"
    }
    else if (x == "count" || x == "Count") {
        x <- "Count"
    }

    Description <- Count <- NULL # to satisfy codetools
    df <- fortify(object, showCategory=showCategory, by=x, ...)

    p <- ggplot(df, aes_string(x = "Description", y = x))
    p <- p + geom_bar(stat = "identity") + coord_flip() + theme_dose(font.size)

    if("pvalue" %in% colnames(p$data)) {
        pvalue <- NULL # to satisfy codetools
        p <- p + aes_string(fill=colorBy) +
            scale_fill_continuous(low="red", high="blue")
    } else {
        p <- p+aes(fill=Description) +
            theme(legend.position="none")
    }
    p <- p + ggtitle(title) + xlab("") + ylab("")
    return(p)
}

