.initial <- function() {
    pos <- 1
    envir <- as.environment(pos)
    assign(".DOSEEnv", new.env(), envir = envir)
    .DOSEEnv <- get(".DOSEEnv", envir = envir)

    tryCatch(utils::data(list="dotbl",
                         package="DOSE"))
    dotbl <- get("dotbl")
    assign("dotbl", dotbl, envir = .DOSEEnv)
    rm(dotbl, envir = .GlobalEnv)

    tryCatch(utils::data(list="mpotbl",
                         package="DOSE"))
    mpotbl <- get("mpotbl")
    assign("mpotbl", mpotbl, envir = .DOSEEnv)
    rm(mpotbl, envir = .GlobalEnv)

    tryCatch(utils::data(list="DOIC",
                         package="DOSE"))
    DOIC <- get("DOIC")
    assign("DOIC", DOIC, envir = .DOSEEnv)
    rm(DOIC, envir = .GlobalEnv)
    
}

check_gene_id <- function(geneList, geneSets) {
    if (all(!names(geneList) %in% unique(unlist(geneSets)))) {
        sg <- unlist(geneSets[1:10])
        sg <- sample(sg, min(length(sg), 6))
        message("--> Expected input gene ID: ", paste0(sg, collapse=','))
        stop("--> No gene can be mapped....")
    }
}


## @importFrom S4Vectors metadata
get_organism <- function(OrgDb) {
    OrgDb <- load_OrgDb(OrgDb)
    ## md <- S4Vectors::metadata(OrgDb)
    ## md[md[,1] == "ORGANISM", 2]
    AnnotationDbi::species(OrgDb)
}


calculate_qvalue <- function(pvals) {
    if (length(pvals) == 0)
        return(numeric(0))

    qobj <- tryCatch(qvalue(pvals, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

    # if (class(qobj) == "qvalue") {
    if (inherits(qobj, "qvalue")) {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    return(qvalues)
}


prepare_relation_df <- function() {
    gtb <- toTable(HDOTERM)
    gtb <- gtb[,1, drop=FALSE]
    gtb <- unique(gtb)

    id <- gtb$do_id
    pid <- mget(id, HDOPARENTS)
    cid <- rep(names(pid), times=sapply(pid, length))

    ptb <- data.frame(id=cid,
                      relationship = 'other',
                      parent = unlist(pid),
                      Ontology = "DO",
                      stringsAsFactors = FALSE)

    dotbl <- merge(gtb, ptb, by.x="doid", by.y="id")
    save(dotbl, file="dotbl.rda", compress="xz")
    invisible(dotbl)
    # mpotbl
    gtb <- toTable(MPO.db::MPOTERM)
    gtb <- gtb[,1, drop=FALSE]
    gtb <- unique(gtb)
    id <- gtb$do_id
    pid <- mget(id, MPO.db::MPOPARENTS)
    cid <- rep(names(pid), times=sapply(pid, length))

    ptb <- data.frame(id=cid,
                      relationship = 'other',
                      parent = unlist(pid),
                      Ontology = "MPO",
                      stringsAsFactors = FALSE)
    
    mpotbl <- merge(gtb, ptb, by.x="mpid", by.y="id")
    save(mpotbl, file="mpotbl.rda", compress="xz")
    invisible(mpotbl)
}




calculate_qvalue <- function(pvals) {
    if (length(pvals) == 0)
        return(numeric(0))

    qobj <- tryCatch(qvalue(pvals, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)
  
    # if (class(qobj) == "qvalue") {
    if (inherits(qobj, "qvalue")) {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    return(qvalues)
}

##' compute information content
##'
##'
##' @title compute information content
##' @param ont one of "DO" and "MPO"
##' @return NULL
##' @importFrom HDO.db HDOTERM
##' @importFrom HDO.db HDOOFFSPRING
##' @importFrom MPO.db MPOMGIDO
##' @importFrom MPO.db MPOANCESTOR
##' @importFrom MPO.db MPOPARENTS
##' @importFrom MPO.db MPOMPMGI
##' @importFrom MPO.db MPOOFFSPRING
##' @importMethodsFrom AnnotationDbi toTable
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
computeIC <- function(ont="DO"){
    if (!exists(".DOSEEnv")) {
        .initial()
    }
    DOSEEnv <- get(".DOSEEnv", envir = .GlobalEnv)
    if (ont == "DO") {
        if (!exists("DO2EG", envir=DOSEEnv)) {
            tryCatch(utils::data(list="DO2EG", package="DOSE"))
            assign("DO2EG", DO2EG, envir = DOSEEnv)
            DO2EG <- get("DO2EG")
            rm(DO2EG, envir = .GlobalEnv)
        }
        DO2EG <- get("DO2EG", envir = DOSEEnv)
        Offsprings <- AnnotationDbi::as.list(HDOOFFSPRING)
    } else {
        eg.do <- toTable(MPOMPMGI)[, c(2,1)]
        colnames(eg.do) <- c("eg", "doid")
        # (2) DOSE:::rebuildAnnoData.internal(eg.do)
        DO2EG <- with(eg.do, split(as.character(eg), as.character(doid)))
        Offsprings <- AnnotationDbi::as.list(MPOOFFSPRING)
    }
    docount <- unlist(lapply(DO2EG, length))
    doids <- names(docount) 
    
    cnt <- docount[doids] + sapply(doids, function(i) sum(docount[Offsprings[[i]]], na.rm=TRUE))
    names(cnt) <- doids
    p <- cnt/sum(docount)

    ## IC of DO terms was quantified as the negative log likelihood.
    IC <- -log(p)
    return(IC)
}


##' provide gene ID, this function will convert to the corresponding DO Terms
##'
##'
##' @title convert Gene ID to DO Terms
##' @param gene entrez gene ID
##' @param organism organism
##' @param ont ont
##' @return DO Terms
##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @export
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
gene2DO <- function(gene, organism = "hsa", ont = "DO") {
    gene <- as.character(gene)
    if (organism == "hsa") {
        if(!exists(".DOSEEnv")) .initial()
        .DOSEEnv <- get(".DOSEEnv", envir=.GlobalEnv)
        if (!exists("EG2DO", envir = .DOSEEnv)) {
            tryCatch(utils::data(list="EG2DO", package="DOSE"))
            EG2DO <- get("EG2DO")
            assign("EG2DO", EG2DO, envir=.DOSEEnv)
            rm(EG2DO, envir=.GlobalEnv)
        }
        EG2DO <- get("EG2DO", envir=.DOSEEnv)
    } else {
        if (ont == "DO") {
            eg.do <- toTable(MPOMGIDO)
            colnames(eg.do) <- c("eg", "doid")
            MPOTERMs <- names(as.list(HDOANCESTOR))             
        } else {
            eg.do <- toTable(MPOMPMGI)[, c(2,1)]
            colnames(eg.do) <- c("eg", "doid")
            MPOTERMs <- names(as.list(MPOANCESTOR))
        }
        EG2DO <- with(eg.do, split(as.character(doid), as.character(eg)))
        EG2DO <- lapply(EG2DO, function(i) unique(i[ i %in% MPOTERMs ]))
        i <- unlist(lapply(EG2DO, function(i) length(i) != 0))
        EG2DO <- EG2DO[i]   
    
    }

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

process_tcss <- getFromNamespace("process_tcss", "GOSemSim")
##' @importClassesFrom GOSemSim GOSemSimDATA
mpodata <- function(processTCSS = FALSE) {
    if (!exists(".DOSEEnv")) .initial()
    DOIC <- new("GOSemSimDATA",
                  ont = "MPO",
                  IC = computeIC(ont = "MPO"))    
    if (processTCSS) {
        message("preparing TCSS data...")
        IC <- DOIC@IC
        DOIC@tcssdata <- process_tcss(ont = "MPO", IC = IC, cutoff = NULL)
    }
    DOIC
}



##' @importClassesFrom GOSemSim GOSemSimDATA
dodata <- function(processTCSS = FALSE) {
    if (!exists(".DOSEEnv")) .initial()
    .DOSEEnv <- get(".DOSEEnv", envir=.GlobalEnv)
    DOIC <- get("DOIC", envir=.DOSEEnv)
    if (processTCSS) {
        message("preparing TCSS data...")
        IC <- DOIC@IC
        DOIC@tcssdata <- process_tcss(ont = "DO", IC = IC, cutoff = NULL)
    }
    DOIC
}

build_dodata <- function() {
    
    DOIC <- new("GOSemSimDATA",
                  ont = "DO",
                  IC = computeIC())
    save(DOIC, file="DOIC.rda", compress="xz")
}


##' rebuilding entrez and DO mapping datasets
##'
##'
##' @title rebuiding annotation data
##' @param file do_rif.human.txt
##' @return NULL
##' @importFrom utils read.delim
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
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

##' @importFrom HDO.db HDOANCESTOR
##' @importFrom HDO.db HDOTERM
##' @importMethodsFrom AnnotationDbi mget
rebuildAnnoData.internal <- function(eg.do) {

    eg <- doid <- NULL # to satisfy codetools

    DO2EG <- with(eg.do, split(as.character(eg), as.character(doid)))
    ## DO2EG <- dlply(eg.do, .(doid), .fun=function(i) i$eg)
    # doids <- toTable(HDOTERM)
    # HDOTERMs <- doids$do_id
    HDOTERMs <- names(as.list(HDOANCESTOR))
    idx <- names(DO2EG) %in% HDOTERMs
    DO2EG <- DO2EG[idx]
    DO2EG <- lapply(DO2EG, function(i) unique(i))
    save(DO2EG, file="DO2EG.rda", compress="xz")

    EG2DO <- with(eg.do, split(as.character(doid), as.character(eg)))
    ## EG2DO <- dlply(eg.do, .(eg), .fun=function(i) i$doid)
    EG2DO <- lapply(EG2DO, function(i) unique(i[ i %in% HDOTERMs ]))

    i <- unlist(lapply(EG2DO, function(i) length(i) != 0))
    EG2DO <- EG2DO[i]
    save(EG2DO, file="EG2DO.rda", compress="xz")

    EG2ALLDO <- lapply(EG2DO,
                       function(i) {
                           ans <- unlist(mget(i, HDOANCESTOR))
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
    ## DO2ALLEG <- dlply(EG2ALLDO.df, .(DO), function(i) as.character(i$EG))
    DO2ALLEG <- with(EG2ALLDO.df, split(as.character(EG), as.character(DO)))
    DO2ALLEG <- lapply(DO2ALLEG, unique)
    save(DO2ALLEG, file="DO2ALLEG.rda", compress="xz")


    ## tryCatch(utils::data(list="DOSEEnv", package="DOSE"))
    ## assign("DO2ALLEG", DO2ALLEG, envir=.DOSEEnv)
    ## assign("EG2ALLDO", EG2ALLDO, envir=.DOSEEnv)
    ## assign("EG2DO", EG2DO, envir=.DOSEEnv)
    ## assign("DO2EG", DO2EG, envir=.DOSEEnv)
    ## save(.DOSEEnv, file="DOSEEnv.rda", compress="xz")
}

## ##' get all entrezgene ID of a specific organism
## ##'
## ##'
## ##' @title getALLEG
## ##' @param organism species
## ##' @return entrez gene ID vector
## ##' @export
## ##' @author Yu Guangchuang
## getALLEG <- function(organism) {
##     annoDb <- getDb(organism)
##     require(annoDb, character.only = TRUE)
##     annoDb <- eval(parse(text=annoDb))
##     eg=keys(annoDb, keytype="ENTREZID")
##     return(eg)
## }


##' mapping gene ID to gene Symbol
##'
##'
##' @title EXTID2NAME
##' @param OrgDb OrgDb
##' @param geneID entrez gene ID
##' @param keytype keytype
##' @return gene symbol
##' @importMethodsFrom AnnotationDbi select
##' @importMethodsFrom AnnotationDbi keys
##' @importMethodsFrom AnnotationDbi columns
##' @importMethodsFrom AnnotationDbi keytypes
##' @importFrom GOSemSim load_OrgDb
##' @export
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
EXTID2NAME <- function(OrgDb, geneID, keytype) {
    OrgDb <- load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
        stop("keytype is not supported...")
    }

    gn.df <- suppressMessages(select(OrgDb, keys=geneID, keytype=keytype, columns="SYMBOL"))
    gn.df <- unique(gn.df)
    colnames(gn.df) <- c("GeneID", "SYMBOL")

    unmap_geneID <- geneID[!geneID %in% gn.df$GeneID]
    if (length(unmap_geneID) != 0) {
        unmap_geneID.df = data.frame(GeneID = unmap_geneID,
                                     SYMBOL = unmap_geneID)
        gn.df <- rbind(gn.df, unmap_geneID.df)
    }

    gn <- gn.df$SYMBOL
    names(gn) <- gn.df$GeneID
    return(gn)
}

## EXTID2NAME <- function(geneID, organism) {
##     if (length(geneID) == 0) {
##         return("")
##     }
##     if (organism == "worm") {
##         organism = "celegans"
##         warning("'worm' is deprecated, please use 'celegans' instead...")
##     }
##     organism <- organismMapper(organism)

##     supported_Org <- getSupported_Org()
##     if (organism %in% supported_Org) {
##         ## kk <- getALLEG(organism)
##         ## unmap_geneID <- geneID[! geneID %in% kk]
##         ## map_geneID <- geneID[geneID %in% kk]

##         ## if (length(map_geneID) == 0) {
##         ##     warning("the input geneID is not entrezgeneID, and cannot be mapped")
##         ##     names(geneID) <- geneID
##         ##     return (geneID)
##         ## }
##         annoDb <- getDb(organism)
##         require(annoDb, character.only = TRUE)
##         annoDb <- eval(parse(text=annoDb))
##         if (organism == "yeast" || organism == "malaria") {
##             gn.df <- select(annoDb, keys=geneID,keytype="ORF", columns="GENENAME")
##         } else if (organism == "arabidopsis") {
##             gn.df <- select(annoDb, keys=geneID,keytype="TAIR", columns="SYMBOL")
##         } else {
##             gn.df <- select(annoDb, keys=geneID,keytype="ENTREZID", columns="SYMBOL")
##         }
##         gn.df <- unique(gn.df)
##         colnames(gn.df) <- c("ENTREZID", "SYMBOL")

##         unmap_geneID <- geneID[!geneID %in% gn.df$ENTREZID]
##         if (length(unmap_geneID) != 0) {
##             unmap_geneID.df = data.frame(ENTREZID= unmap_geneID, SYMBOL=unmap_geneID)
##             gn.df <- rbind(gn.df, unmap_geneID.df)
##         }

##         gn <- gn.df$SYMBOL
##         names(gn) <- gn.df$ENTREZID
##         ##gn <- unique(gn[!is.na(gn)])
##     } else {
##         oldwd <- getwd()
##         if(organism == "D39") {
##             dir <- system.file("extdata/D39/", package="clusterProfiler")
##             setwd(dir)
##         }
##         if(organism == "M5005") {
##             dir <- system.file("extdata/M5005/", package="clusterProfiler")
##             setwd(dir)
##         }

##         if (file.exists("geneTable.rda")) {
##             geneTable <- NULL # to satisfy codetools
##             load("geneTable.rda")
##             idx <- geneTable$GeneID %in% geneID
##             eg.gn <- geneTable[idx, c("GeneID", "GeneName", "Locus")]
##             eg.gn[eg.gn[,2] == "-",2] <- eg.gn[eg.gn[,2] == "-",3]
##             ##eg.gn <- eg.gn[,c(1,2)]
##             gn <- eg.gn$GeneName
##             names(gn) <- as.character(eg.gn$GeneID)
##             setwd(oldwd)
##         } else {
##             setwd(oldwd)
##             warning("Have no annotation found for the input geneID")
##             return(geneID)
##         }
##     }
##     return(gn)
## }



is.sorted <- function(x, decreasing=TRUE) {
    all( sort(x, decreasing=decreasing) == x )
}

getGeneSet <- function(USER_DATA) {
    if (inherits(USER_DATA, "environment")) { 
        res <- get("PATHID2EXTID", envir = USER_DATA)
    } else if (inherits(USER_DATA, "GSON")) {
        gsid2gene <- USER_DATA@gsid2gene
        res <- split(gsid2gene$gene, gsid2gene$gsid) 
    } else {
        stop("not supported")
    }
    return(res)
}


##' @importFrom ggplot2 facet_grid
##' @export
ggplot2::facet_grid
