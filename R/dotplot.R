##' dot plot of enrichResult
##'
##' 
##' @title dotplot
##' @param object an instance of enrichResult
##' @param x variable for x axis
##' @param colorBy one of 'pvalue', 'p.adjust' and 'qvalue'
##' @param showCategory number of category
##' @param font.size font size
##' @param title plot title
##' @return ggplot object
##' @importFrom ggplot2 fortify
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 scale_color_gradient
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @export
##' @author Guangchuang Yu
dotplot <- function(object, x="geneRatio", colorBy="p.adjust", showCategory=10, font.size=12, title="") {
    if (! is(object, "enrichResult")) {
        stop("object should be an instance of 'enrichResult'...")
    }
    colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue"))
    if (x == "geneRatio" || x == "GeneRatio") {
        x <- "GeneRatio"
        size <- "Count"
    } else if (x == "count" || x == "Count") {
        x <- "Count"
        size <- "GeneRatio"
    } else {
        stop("x should be geneRatio or count...")
    }
    df <- fortify(object, showCategory = showCategory)
    gsize <- as.numeric(sub("/\\d+$", "", as.character(df$GeneRatio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(df$GeneRatio)))
    df$GeneRatio <- gsize/gcsize

    idx <- order(df$GeneRatio, decreasing = FALSE)
    df$Description <- factor(df$Description, level=df$Description[idx])
    ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
        geom_point() + scale_color_gradient(low="red", high="blue") +
            ylab("") + ggtitle(title) + theme_dose(font.size)
}
