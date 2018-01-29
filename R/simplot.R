##' plotting similarity matrix
##'
##'
##' @title simplot
##' @param sim similarity matrix
##' @param xlab xlab
##' @param ylab ylab
##' @param color.low color of low value
##' @param color.high color of high value
##' @param labs logical, add text label or not
##' @param digits round digit numbers
##' @param labs.size lable size
##' @param font.size font size
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 geom_text
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
simplot <- function(sim, xlab="", ylab="", color.low="white", color.high="red", labs=TRUE, digits=2, labs.size=3, font.size=14){
    sim.df <- as.data.frame(sim)
    ## if(readable == TRUE) {
    ##     rownames(sim.df) <- TERM2NAME(rownames(sim.df))
    ##     colnames(sim.df) <- TERM2NAME(colnames(sim.df))
    ## }
    rn <- row.names(sim.df)

    sim.df <- cbind(ID=rownames(sim.df), sim.df)
    sim.df <- melt(sim.df)

    sim.df[,1] <- factor(sim.df[,1], levels=rev(rn))
    if (labs == TRUE) {
        ## lbs <- c(apply(round(sim, digits), 2, as.character))
        sim.df$label <- as.character(round(sim.df$value, digits))
    }
    variable <- ID <- value <- label <- NULL ## to satisfy codetools
    if (labs == TRUE)
        p <- ggplot(sim.df, aes(variable, ID, fill=value, label=label))
    else
        p <- ggplot(sim.df, aes(variable, ID, fill=value))

     p <- p + geom_tile(color="black")+
            scale_fill_gradient(low=color.low, high=color.high) +
                scale_x_discrete(expand=c(0,0)) +
                    scale_y_discrete(expand=c(0,0))+
                        theme(axis.ticks=element_blank())
    if (labs == TRUE)
        p <- p+geom_text(size=labs.size)
    p <- p+theme_dose(font.size)
    p <- p + theme(axis.text.x=element_text(hjust=0, angle=-90)) +
        theme(axis.text.y=element_text(hjust=0))
    p <- p+theme(legend.title=element_blank())
    ##geom_point(aes(size=value))
    p <- p+xlab(xlab)+ylab(ylab)

    ## if (readable == TRUE) {
    ##     p <- p + theme(axis.text.y = element_text(hjust=1))
    ## }
    p <- p + theme(axis.text.x = element_text(vjust=0.5))
    return(p)
}


##' ggplot theme of DOSE
##'
##' @title theme_dose
##' @param font.size font size
##' @return ggplot theme
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @examples
##' library(ggplot2)
##' qplot(1:10) + theme_dose()
##' @export
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
