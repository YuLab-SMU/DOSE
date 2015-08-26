.initial <- function() {
    assign("DOSEEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("ICEnv", new.env(), .GlobalEnv)

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
    EG2DO <- get("EG2DO", envir=DOSEEnv)
    DO <- EG2DO[[gene]]
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
##' @author Guangchuang Yu \url{http://ygc.name}
rebuildAnnoData <- function(file) {
    ##
    ## do_rif.human.txt was downloaded from
    ## http://projects.bioinformatics.northwestern.edu/do_rif/
    ##
    
    ## do.rif <- read.delim2(file, sep="\t", stringsAsFactors=F, header=F)
    ## eg.do <- do.rif[,c(1,5)]

    ## new file
    ## IDMappings.txt from
    ## http://doa.nubic.northwestern.edu/pages/download.php
    domapping <- read.delim(file, stringsAsFactors=F)
    eg.do <- domapping[,c(2,1)]
    colnames(eg.do) <- c("eg", "doid")
    eg.do$doid <- paste("DOID:", eg.do$doid, sep="")
    
    rebuildAnnoData.internal(eg.do)
}

##' @importFrom plyr dlply
##' @importFrom plyr .
##' @importFrom DO.db DOANCESTOR
##' @importFrom DO.db DOTERM
##' @importMethodsFrom AnnotationDbi mget
rebuildAnnoData.internal <- function(eg.do) {
    
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


    tryCatch(utils::data(list="DOSEEnv", package="DOSE"))    
    assign("DO2ALLEG", DO2ALLEG, envir=DOSEEnv)
    assign("EG2ALLDO", EG2ALLDO, envir=DOSEEnv)
    assign("EG2DO", EG2DO, envir=DOSEEnv)
    assign("DO2EG", DO2EG, envir=DOSEEnv)
    save(DOSEEnv, file="DOSEEnv.rda", compress="xz")
}

##' get all entrezgene ID of a specific organism
##'
##'
##' @title getALLEG
##' @param organism species
##' @return entrez gene ID vector
##' @export
##' @author Yu Guangchuang
getALLEG <- function(organism) {
    annoDb <- getDb(organism)
    require(annoDb, character.only = TRUE)
    annoDb <- eval(parse(text=annoDb))
    eg=keys(annoDb, keytype="ENTREZID")
    return(eg)
}


##' mapping gene ID to gene Symbol
##'
##'
##' @title EXTID2NAME
##' @param geneID entrez gene ID
##' @param organism one of "human", "mouse" and "yeast"
##' @return gene symbol
##' @importMethodsFrom AnnotationDbi select
##' @importMethodsFrom AnnotationDbi keys
##' @importMethodsFrom AnnotationDbi columns
##' @importFrom GOSemSim getSupported_Org
##' @importFrom GOSemSim getDb
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
EXTID2NAME <- function(geneID, organism) {
    if (length(geneID) == 0) {
        return("")
    }
    if (organism == "worm") {
        organism = "celegans"
        warning("'worm' is deprecated, please use 'celegans' instead...")
    }
    organism <- organismMapper(organism)
    
    supported_Org <- getSupported_Org()
    if (organism %in% supported_Org) {
        ## kk <- getALLEG(organism)
        ## unmap_geneID <- geneID[! geneID %in% kk]
        ## map_geneID <- geneID[geneID %in% kk]

        ## if (length(map_geneID) == 0) {
        ##     warning("the input geneID is not entrezgeneID, and cannot be mapped")
        ##     names(geneID) <- geneID
        ##     return (geneID)
        ## }
        annoDb <- getDb(organism)
        require(annoDb, character.only = TRUE)
        annoDb <- eval(parse(text=annoDb))
        if (organism == "yeast" || organism == "malaria") {
            gn.df <- select(annoDb, keys=geneID,keytype="ORF", columns="GENENAME")
        } else if (organism == "arabidopsis") {
            gn.df <- select(annoDb, keys=geneID,keytype="TAIR", columns="SYMBOL")
        } else {
            gn.df <- select(annoDb, keys=geneID,keytype="ENTREZID", columns="SYMBOL")
        }
        gn.df <- unique(gn.df)
        colnames(gn.df) <- c("ENTREZID", "SYMBOL")

        unmap_geneID <- geneID[!geneID %in% gn.df$ENTREZID]
        if (length(unmap_geneID) != 0) {
            unmap_geneID.df = data.frame(ENTREZID= unmap_geneID, SYMBOL=unmap_geneID)
            gn.df <- rbind(gn.df, unmap_geneID.df)
        }

        gn <- gn.df$SYMBOL
        names(gn) <- gn.df$ENTREZID
        ##gn <- unique(gn[!is.na(gn)])
    } else {
        oldwd <- getwd()
        if(organism == "D39") {
            dir <- system.file("extdata/D39/", package="clusterProfiler")
            setwd(dir)
        }
        if(organism == "M5005") {
            dir <- system.file("extdata/M5005/", package="clusterProfiler")
            setwd(dir)
        }
        
        if (file.exists("geneTable.rda")) {
            geneTable <- NULL # to satisfy codetools
            load("geneTable.rda")
            idx <- geneTable$GeneID %in% geneID
            eg.gn <- geneTable[idx, c("GeneID", "GeneName", "Locus")]
            eg.gn[eg.gn[,2] == "-",2] <- eg.gn[eg.gn[,2] == "-",3]
            ##eg.gn <- eg.gn[,c(1,2)]
            gn <- eg.gn$GeneName
            names(gn) <- as.character(eg.gn$GeneID)
            setwd(oldwd)
        } else {
            setwd(oldwd)
            warning("Have no annotation found for the input geneID")
            return(geneID)
        }
    }
    return(gn)
}



organismMapper <- function(organism) {
    ## to satisfy the change of enrichKEGG
    
    if (organism == "aga") {
        species <- "anopheles"
    } else if (organism == "ath") {
        species <- "arabidopsis"
    } else if (organism == "bta") {
        species <- "bovine"
    } else if (organism == "cfa") {
        species <- "canine"
    } else if (organism == "gga") {
        species <- "chicken"
    } else if (organism == "ptr") {
        species <- "chipm"
    } else if (organism == "eco") {
        species <- "ecolik12"
    } else if (organism == "ecs") {
        species <- "ecsakai"
    } else if (organism == "dme") {
        species <- "fly"
    } else if (organism == "hsa") {
        species <- "human"
    } else if (organism == "pfa") {
        species <- "malaria"
    } else if (organism == "mmu") {
        species <- "mouse"
    } else if (organism == "ssc") {
        species <- "pig"
    } else if (organism == "rno") {
        species <- "rat"
    } else if (organism == "mcc") {
        species <- "rhesus"
    } else if (organism == "cel") {
        species <- "worm"
    } else if (organism == "xla") {
        species <- "xenopus"
    } else if (organism == "sce") {
        species <- "yeast"
    } else if (organism == "dre") {
        species <- "zebrafish"
    } else {
        species <- organism
    }
    return(species)
}
