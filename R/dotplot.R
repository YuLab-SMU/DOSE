##' @importFrom ggplot2 fortify
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 scale_color_gradient
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
dotplot_internal <- function(object, x="geneRatio", colorBy="p.adjust", showCategory=10, split=NULL, font.size=12, title="") {
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
    df <- fortify(object, showCategory = showCategory, split=split)
    ## already parsed in fortify
    ## df$GeneRatio <- parse_ratio(df$GeneRatio)

    idx <- order(df$GeneRatio, decreasing = FALSE)
    df$Description <- factor(df$Description, levels=df$Description[idx])
    ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
        geom_point() + scale_color_gradient(low="red", high="blue") +
            ylab("") + ggtitle(title) + theme_dose(font.size)
}

