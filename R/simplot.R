##' plotting similarity matrix
##'
##'
##' @title simplot
##' @param sim similarity matrix
##' @param xlab xlab
##' @param ylab ylab
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_gradient
##' @importFrom ggplot2 scale_x_discrete
##' @importFrom ggplot2 scale_y_discrete
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom reshape2 melt
##' @export
##' @author Yu Guangchuang
simplot <- function(sim, xlab="", ylab="") {
    sim.df <- as.data.frame(sim)
    sim.df <- cbind(ID=rownames(sim.df), sim.df)
    sim.df <- melt(sim.df)

    sim.df[,1] <- factor(sim.df[,1], levels=rev(rownames(sim)))

    variable <- ID <- value <- NULL ## to satisfy codetools
    p <- ggplot(sim.df, aes(variable, ID, fill=value)) +
	geom_tile(color="black")+
            scale_fill_gradient(low="white", high="steelblue") +
                scale_x_discrete(expand=c(0,0)) +
                    scale_y_discrete(expand=c(0,0))+
                        theme(axis.ticks=element_blank())

    p <- p + theme(axis.text.x=element_text(size=12, hjust=0, angle=-90)) +
	theme(axis.text.y=element_text(size=12, hjust=0))
    p <- p+theme(legend.title=element_blank())
    ##geom_point(aes(size=value))
    p <- p+xlab(xlab)+ylab(ylab)
    return(p)
}
