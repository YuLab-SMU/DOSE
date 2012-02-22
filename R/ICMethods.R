##' compute information content
##'
##'
##' @title compute information content
##' @param ont "DO"
##' @param organism "human"
##' @return NULL
##' @importFrom DO.db DOTERM
##' @importFrom DO.db DOOFFSPRING
##' @importMethodsFrom AnnotationDbi toTable
##' @author Guangchuang Yu \url{http://ygc.name}
computeIC <- function(ont="DO", organism="human"){
    doids <- toTable(DOTERM)
    doterms <- doids$do_id
    docount <- table(doterms)
    doids <- names(docount)  #unique(doterms)
    cnt <- sapply(doids,function(x){
        n=docount[get(x, DOOFFSPRING)]
        docount[x]+sum(n[!is.na(n)])
    })
    names(cnt) <- doids
    p <- cnt/sum(docount)

    ## IC of DO terms was quantified as the negative log likelihood.
    IC <- -log(p)
    fname <- paste(paste("Info_Contents",
                         organism,
                         ont,
                         sep="_"),
                   ".rda",
                   sep="")
    save(IC, file=fname)
}

##' Load Information Content data to DOSEEnv environment
##'
##'
##' @title Load IC data
##' @param organism "human"
##' @param ont "DO"
##' @return NULL
##' @author Guangchuang Yu \url{http://ygc.name}
loadICdata <- function(organism, ont) {
    if(!exists("ICEnv")) .initial()
    fname <- paste("Info_Contents",
                   organism,
                   ont,
                   sep="_")
    tryCatch(utils::data(list=fname,
                         package="DOSE"))
    IC <- get("IC")
    org.ont.IC <- paste(organism,
                        ont,
                        "IC",
                        sep="")
    assign(eval(org.ont.IC),
           IC,
           envir=ICEnv)
    rm (IC)
}


##' Information Content Based Methods for semantic similarity measuring
##'
##' implemented for methods proposed by Resnik, Jiang, Lin and Schlicker.
##' @title information content based methods
##' @param ID1 Ontology Term
##' @param ID2 Ontology Term
##' @param ont Ontology
##' @param method one of "Resnik", "Jiang", "Lin" and "Rel".
##' @param organism one of supported species
##' @return semantic similarity score
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
infoContentMethod <- function(ID1,
                              ID2,
                              ont="DO",
                              method,
                              organism="human") {
    if(!exists("ICEnv")) {
        .initial()
    }

    org.ont.IC <- paste(organism,
                        ont,
                        "IC",
                        sep="")

    if(!exists(org.ont.IC, envir=ICEnv)) {
        loadICdata(organism, ont)
    }
    IC <- get(org.ont.IC, envir=ICEnv)

    ## more specific term, larger IC value.
    ## Normalized, all divide the most informative IC.
    ## all IC values range from 0(root node) to 1(most specific node)
    mic <- max(IC[IC!=Inf])

    if (ont == "DO") {
        topNode <- "DOID:4"
    } else {
        topNode <- "all"
    }

    IC[topNode] = 0

    ic1 <- IC[ID1]/mic
    ic2 <- IC[ID2]/mic

    if (ic1 == 0 || ic2 == 0)
        return (NA)

    ONTANCESTOR <- .getAncestors(ont)
    ancestor1 <- get(ID1, ONTANCESTOR)
    ancestor2 <- get(ID2, ONTANCESTOR)
    if (ID1 == ID2) {
        commonAncestor <- ID1
    } else if (ID1 %in% ancestor2) {
        commonAncestor <- ID1
    } else if (ID2 %in% ancestor1) {
        commonAncestor <- ID2
    } else {
        commonAncestor <- intersect(ancestor1, ancestor2)
    }
    if (length(commonAncestor) == 0) return (NA)

    ##Information Content of the most informative common ancestor (MICA)
    mica <- max(IC[commonAncestor])/mic

    ## IC is biased
    ## because the IC of a term is dependent of its children but not on its parents.
    sim <- switch(method,
                  Resnik = mica, ## Resnik does not consider how distant the terms are from their common ancestor.
                  ## Lin and Jiang take that distance into account.
                  Lin = 2*mica/(ic1+ic2),
                  Jiang = 1 - min(1, -2*mica + ic1 + ic2),
                  Rel = 2*mica/(ic1+ic2)*(1-exp(-mica*mic))  ## mica*mic equals to the original IC value. and exp(-mica*mic) equals to the probability of the term's occurence.
                  )
    return (sim)
}


##' @importFrom DO.db DOANCESTOR
.getAncestors <- function(ont) {
    Ancestors <- switch(ont,
                        DO = DOANCESTOR
                        )
    return(Ancestors)
}
