
##' Method Wang for semantic similarity measuring
##'
##'
##' @title wangMethod
##' @param ID1 Ontology Term
##' @param ID2 Ontology Term
##' @param ont Ontology
##' @param weight.isa weight of isa relationship
##' @param weight.partof weight of partof relationship
##' @param weight.do weight of DO
##' @return semantic similarity score
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
wangMethod <- function(ID1,
                       ID2,
                       ont="DO",
                       weight.isa=0.8,
                       weight.partof=0.6,
                       weight.do=0.7) {

    if (ID1 == ID2)
        return (sim=1)

    sv.a <- 1
    sv.b <- 1
    sw <- 1
    names(sv.a) <- ID1
    names(sv.b) <- ID2

    Parents <- .getParents(ont)
    sv.a <- .SemVal(ID1,
                    ont,
                    Parents,
                    sv.a,
                    sw,
                    weight.isa,
                    weight.partof,
                    weight.do)

    sv.b <- .SemVal(ID2,
                    ont,
                    Parents,
                    sv.b,
                    sw,
                    weight.isa,
                    weight.partof,
                    weight.do)

    sv.a <- .uniqsv(sv.a)
    sv.b <- .uniqsv(sv.b)

    idx <- intersect(names(sv.a), names(sv.b))
    inter.sva <- unlist(sv.a[idx])
    inter.svb <- unlist(sv.b[idx])
    if (is.null(inter.sva) ||
        is.null(inter.svb) ||
        length(inter.sva) == 0 ||
        length(inter.svb) ==0) {
        sim <- NA
    } else {
        sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
    }
    return(sim)
}

.uniqsv <- function(sv) {
    sv <- unlist(sv)
    una <- unique(names(sv))
    sv <- unlist(sapply(una,
                        function(x) {
                            max(sv[names(sv)==x])
                        }
                        )
                 )
    return (sv)
}


##' @importFrom DO.db DOPARENTS
.getParents <- function(ont) {
	Parents <- switch(ont,
		DO = DOPARENTS
	)
	return(Parents)
}

.SemVal_internal <- function(ID,
                             ont,
                             Parents,
                             sv,
                             w,
                             weight.isa,
                             weight.partof,
                             weight.do) {

    if (!exists(ID, Parents)) {
        return(NA)
    }
    p <- get(ID, Parents)
    ##p <- unlist(p[[1]])
    if (length(p) == 0 || is.na(p)) {
        ##warning(ID, " may not belong to Ontology ", ont)
        return(NA)
    }

    old.w <- w
    if (ont == "DO") {
        topNode <- "DOID:4"
    } else {
        relations <- names(p)
        topNode <- "all"
    }

    for (i in 1:length(p)) {
        if (ont == "DO") {
            w <- old.w * weight.do
        } else {
            if (relations[i] == "is_a") {
                w <- old.w * weight.isa
            } else {
                w <- old.w * weight.partof
            }
        }
        names(w) <- p[i]
        sv <- c(sv,w)
        if (p[i] != topNode) {
            sv <- .SemVal_internal(p[i], ont, Parents, sv, w, weight.isa, weight.partof, weight.do)
        }
    }
    return (sv)
}

.SemVal <- function(ID,
                    ont,
                    Parents,
                    sv,
                    w,
                    weight.isa,
                    weight.partof,
                    weight.do) {
    ##	if(!exists("SemSimCache")) return(.SemVal_internal(ID, ont, Parents, sv, w, weight.isa, weight.partof, weight.do))
    if(!exists("SemSimCache")) {
        .initial()
    }
    ID.ont <- paste(ID, ont, sep=".")
    if (!exists(ID.ont, envir=SemSimCache)) {
        value <- .SemVal_internal(ID,
                                  ont,
                                  Parents,
                                  sv,
                                  w,
                                  weight.isa,
                                  weight.partof,
                                  weight.do)

        assign(ID.ont,
               value,
               envir=SemSimCache)
                                        #cat("recompute ", ID, value, "\n")
    }
    else{
                                        #cat("cache ", ID, get(ID, envir=SemSimCache), "\n")
    }
    return(get(ID.ont,
               envir=SemSimCache)
           )
}
