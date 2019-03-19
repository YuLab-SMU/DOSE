##' parse character ratio to double value, such as 1/5 to 0.2
##'
##' 
##' @title parse_ratio
##' @param ratio character vector of ratio to parse
##' @return A numeric vector (double) of parsed ratio
##' @export
##' @author Guangchuang Yu
parse_ratio <- function(ratio) {
    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator/denominator)
}

