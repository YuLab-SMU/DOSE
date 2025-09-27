get_dose_env <- function() {
    if (!exists(".DOSEEnv")) {
        .initial()
    }
    get(".DOSEEnv")
}

.initial <- function() {
    pos <- 1
    envir <- as.environment(pos)
    assign(".DOSEEnv", new.env(), envir = envir) 
}

# https://github.com/YuLab-SMU/ReactomePA/issues/43
check_gene_id <- function(geneList, geneSets) {
    if (all(!names(geneList) %in% unique(unlist(geneSets)))) {
        sg <- unlist(geneSets[1:10])
        sg <- sample(sg, min(length(sg), 6))
        message("--> Expected input gene ID: ", paste0(sg, collapse=','))
        message("--> No gene can be mapped....")
        return(FALSE)
    }
    return(TRUE)
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
##' @param ont one of "DO", "HPO" and "MPO"
##' @return NULL
##' @importMethodsFrom AnnotationDbi toTable
##' @author Guangchuang Yu \url{https://yulab-smu.top}
computeIC <- function(ont="HDO"){
    DO2EG <- get_ont2gene(ont)
    Offsprings <- GOSemSim:::getOffsprings(ont)
    
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
##' @author Guangchuang Yu \url{https://yulab-smu.top}
gene2DO <- function(gene, organism = "hsa", ont = "HDO") {
    gene <- as.character(gene)

    EG2DO <- get_gene2ont(ont)

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
semdata <- function(processTCSS = FALSE, ont = "HDO") {
    IC <- new("GOSemSimDATA",
                ont = ont,
                IC = computeIC(ont = ont))

    if (processTCSS) {
        IC <- IC@IC
        IC@tcssdata <- process_tcss(ont = ont, IC = IC, cutoff = NULL)
    }

    IC
}

semdata2 <- memoise::memoise(semdata)


get_ont2gene <- function(ontology, output = "list") {
    gene2ont <- get_gene2ont(ontology, output = "data.frame")
    if (output == "data.frame") {
        return(gene2ont[, 2:1])
    }

    split(as.character(gene2ont[,1]), as.character(gene2ont[,2]))
}

get_gene2ont <- function(ontology, output = "list") {
    ont2gene <- GOSemSim:::get_onto_data(ontology, table = "ont2gene", output = 'data.frame')
    anc <- GOSemSim:::getAncestors(ontology)
    idx <- ont2gene[,1] %in% names(anc)
    ont2gene <- unique(ont2gene[idx, ])

    if (output == "data.frame") {
        return(ont2gene[, 2:1])
    }

    split(as.character(ont2gene[,1]), as.character(ont2gene[,2]))
}

get_gene2allont <- function(ontology, output = "list") {
    GOSemSim:::get_onto_data(ontology, table = "gene2allont", output = output)
}

get_ont2allgene <- function(ontology, output = "list") {
    gene2allont <- GOSemSim:::get_onto_data(ontology, table = "gene2allont", output = "data.frame")
    if (output == "data.frame") {
        return(gene2allont[, 2:1])
    }

    split(as.character(gene2allont[,1]), as.character(gene2allont[,2]))
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
##' @author Guangchuang Yu \url{https://yulab-smu.top}
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
