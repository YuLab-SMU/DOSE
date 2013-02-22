.initial <- function() {
    assign("DOSEEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("ICEnv", new.env(), .GlobalEnv)

    ## tryCatch(utils::data(list="DO2ALLEG", package="DOSE"))
    ## assign("DO2ALLEG", DO2ALLEG, envir=DOSEEnv)

    ## tryCatch(utils::data(list="EG2ALLDO", package="DOSE"))
    ## assign("EG2ALLDO", EG2ALLDO, envir=DOSEEnv)

    ## tryCatch(utils::data(list="EG2DOLite", package="DOSE"))
    ## assign("EG2DOLite", EG2DOLite, envir=DOSEEnv)

    ## tryCatch(utils::data(list="DOLite2EG", package="DOSE"))
    ## assign("DOLite2EG", DOLite2EG, envir=DOSEEnv)

    ## tryCatch(utils::data(list="DOLiteTerm", package="DOSE"))
    ## assign("DOLiteTerm", DOLiteTerm, envir=DOSEEnv)

    tryCatch(utils::data(list="DOSEEnv", package="DOSE"))
}


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


##' provide gene ID, this function will convert to the corresponding DO Terms
##'
##'
##' @title convert Gene ID to DO Terms
##' @param gene entrez gene ID
##' @return DO Terms
##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @author Guangchuang Yu \url{http://ygc.name}
gene2DO <- function(gene) {
    if(!exists("DOSEEnv")) .initial()
    EG2ALLDO <- get("EG2ALLDO", envir=DOSEEnv)
    DO <- EG2ALLDO[[gene]]
    DO <- unlist(DO)
    if (is.null(DO)) {
        return(NA)
    }
    if (sum(!is.na(DO)) == 0) {
        return(NA)
    }
    DO <- DO[!is.na(DO)]
    if (length(DO) == 0) {
        return(NA)
    }
    return(DO)
}

##' rebuilding entrez and DO mapping datasets
##'
##'
##' @title rebuiding annotation data
##' @param file do_rif.human.txt
##' @return NULL
##' @importFrom plyr dlply
##' @importFrom plyr .
##' @importFrom DO.db DOANCESTOR
##' @importFrom DO.db DOTERM
##' @importMethodsFrom AnnotationDbi mget
##' @author Guangchuang Yu \url{http://ygc.name}
rebuildAnnoData <- function(file) {
                                        #
                                        # do_rif.human.txt was downloaded from
                                        # http://projects.bioinformatics.northwestern.edu/do_rif/
                                        #
    do.rif <- read.delim2(file, sep="\t", stringsAsFactors=F, header=F)
    eg.do <- do.rif[,c(1,5)]
    colnames(eg.do) <- c("eg", "doid")

    eg <- doid <- NULL # to satisfy codetools

    DO2EG <- dlply(eg.do, .(doid), .fun=function(i) i$eg)
    doids <- toTable(DOTERM)
    doterms <- doids$do_id
    idx <- names(DO2EG) %in% doterms
    DO2EG <- DO2EG[idx]
    DO2EG <- lapply(DO2EG, function(i) unique(i))
    save(DO2EG, file="DO2EG.rda", compress="xz")


    EG2DO <- dlply(eg.do, .(eg), .fun=function(i) i$doid)
    EG2DO <- lapply(EG2DO, function(i) unique(i[ i %in% doterms ]))

    i <- unlist(lapply(EG2DO, function(i) length(i) != 0))
    EG2DO <- EG2DO[i]
    save(EG2DO, file="EG2DO.rda", compress="xz")

    EG2ALLDO <- lapply(EG2DO,
                       function(i) {
                           ans <- unlist(mget(i, DOANCESTOR))
                           ans <- ans[ !is.na(ans) ]
                           ans <- c(i, ans)
                           ans <- unique(ans)
                           return(ans)
                       })
    save(EG2ALLDO, file="EG2ALLDO.rda", compress="xz")

    len <- lapply(EG2ALLDO,length)
    EG2ALLDO.df <- data.frame(EG=rep(names(EG2ALLDO), times=len),
                              DO=unlist(EG2ALLDO))
    DO <- NULL ## satisfy code tools
    DO2ALLEG <- dlply(EG2ALLDO.df, .(DO), function(i) as.character(i$EG))
    DO2ALLEG <- lapply(DO2ALLEG, unique)
    save(DO2ALLEG, file="DO2ALLEG.rda", compress="xz")
}

##' mapping gene ID to gene Symbol
##'
##'
##' @title EXTID2NAME
##' @param geneID entrez gene ID
##' @param organism one of "human", "mouse" and "yeast"
##' @return gene symbol
##' @importMethodsFrom AnnotationDbi select
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
EXTID2NAME <- function(geneID, organism) {
    if (length(geneID) == 0) {
        return("")
    }

    supported_Org <- c("human", "mouse", "yeast", "zebrafish", "celegans")
    if (organism %in% supported_Org) {
        annoDb <- switch(organism,
                         human = "org.Hs.eg.db",
                         mouse = "org.Mm.eg.db",
                         yeast = "org.Sc.sgd.db",
                         zebrafish = "org.Dr.eg.db",
                         celegans="org.Ce.eg.db"
                         )
        require(annoDb, character.only=TRUE)
        annoDb <- eval(parse(text=annoDb))
        geneID <- as.character(geneID)

        kk=keys(annoDb, keytype="ENTREZID")
        unmap_geneID <- geneID[! geneID %in% kk]
        map_geneID <- geneID[geneID %in% kk]

        if (length(map_geneID) == 0) {
            warning("the input geneID is not entrezgeneID, and cannot be mapped")
            names(geneID) <- geneID
            return (geneID)
        }
        gn.df <- select(annoDb, keys=geneID,cols="SYMBOL")
        gn.df <- unique(gn.df)

        if (length(unmap_geneID) != 0) {
            unmap_geneID.df = data.frame(ENTREZID= unmap_geneID, SYMBOL=unmap_geneID)
            gn.df <- rbind(gn.df, unmap_geneID.df)
        }

        gn <- gn.df$SYMBOL
        names(gn) <- gn.df$ENTREZID
        ##gn <- unique(gn[!is.na(gn)])
    } else {
        if (file.exists("geneTable.rda")) {
            geneTable <- NULL # to satisfy codetools
            load("geneTable.rda")
            idx <- geneTable$GeneID %in% geneID
            eg.gn <- geneTable[idx, c("GeneID", "GeneName", "Locus")]
            eg.gn[eg.gn[,2] == "-",2] <- eg.gn[eg.gn[,2] == "-",3]
            ##eg.gn <- eg.gn[,c(1,2)]
            gn <- eg.gn$GeneName
            names(gn) <- as.character(eg.gn$GeneID)
        } else {
            warning("Have no annotation found for the input geneID")
            return(geneID)
        }
    }
    return(gn)
}
