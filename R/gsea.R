##' generic function for gene set enrichment analysis
##'
##'
##' @title gsea
##' @param geneList order ranked geneList
##' @param geneSets gene sets
##' @param setType Type of geneSet
##' @param organism organism
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param pvalueCutoff p value Cutoff
##' @param pAdjustMethod  p value adjustment method
##' @param verbose print message or not
##' @param ... additional parameter
##' @return gseaResult object
##' @importFrom plyr ldply
##' @importFrom parallel detectCores
##' @importFrom parallel mclapply
##' @export
##' @author Yu Guangchuang
gsea <- function(geneList,
                 geneSets,
                 setType,
                 organism,
                 exponent,
                 nPerm,
                 minGSSize,
                 pvalueCutoff,
                 pAdjustMethod,
                 verbose, ...) {

    class(setType) <- "character"

    ## index of geneSets in used.
    ## logical
    geneSets <- sapply(geneSets, intersect, names(geneList))
    gs.idx <- sapply(geneSets, length) > minGSSize
    nGeneSet <- sum(gs.idx)

    if ( nGeneSet == 0 )
        return (NULL)

    selected.gs <- geneSets[gs.idx]
    if (verbose)
        print("calculating observed enrichment scores...")
    observedScore <- sapply(selected.gs, function(gs)
                            gseaScores(geneSet=gs,
                                       geneList=geneList,
                                       exponent=exponent)
                            )

    if(verbose) {
        print("calculating permutation scores...")
        pb <- txtProgressBar(min=0, max=nGeneSet, style=3)
    }
    if(Sys.info()[1] == "Windows") {
        permScores <- t(sapply(seq_along(selected.gs), function(i) {
            if(verbose)
                setTxtProgressBar(pb, i)
            perm.gseaEScore(geneList=geneList,
                            geneSet=selected.gs[[i]],
                            nPerm=nPerm,
                            exponent=exponent)
        }))
    } else {
        permScores <- mclapply(seq_along(selected.gs), function(i) {
            if(verbose)
                setTxtProgressBar(pb, i)
            perm.gseaEScore(geneList=geneList,
                            geneSet=selected.gs[[i]],
                            nPerm=nPerm,
                            exponent=exponent)
        },
                               mc.cores=detectCores()
                               )
        permScores <- ldply(permScores)
        permScores <- as.matrix(permScores)
    }

    if(verbose)
        close(pb)

    rownames(permScores) <- names(selected.gs)

    if (verbose)
        print("calculating p values...")
    pvals <- sapply(seq_along(observedScore), function(i) {
        if( is.na(observedScore[i]) ) {
            NA
        } else if ( observedScore[i] == 0 ) {
            1
        } else if ( observedScore[i] > 0 ) {
            (sum(permScores[i, ] >= observedScore[i]) +1) / (nPerm+1)
        } else { # observedScore[i] < 0
            (sum(permScores[i, ] <= observedScore[i]) +1) / (nPerm+1)
        }
    })
    p.adj <- p.adjust(pvals, method=pAdjustMethod)
    qobj <- qvalue(pvals, lambda=0.05, pi0.method="bootstrap")
    if (class(qobj) == "qvalue") {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }

    gs.name <- names(selected.gs)
    class(gs.name) <- setType
    Description <- TERM2NAME(gs.name, organism, ...)

    params <- list(setType = setType,
                   organism = organism,
                   pvalueCutoff = pvalueCutoff,
                   nPerm = nPerm,
                   pAdjustMethod = pAdjustMethod,
                   exponent = exponent,
                   minGSSize = minGSSize
                   )

    res <- data.frame(
                      ID = as.character(gs.name),
                      Description = Description,
                      setSize = sapply(selected.gs, length),
                      enrichmentScore = observedScore,
                      pvalue = pvals,
                      p.adjust = p.adj,
                      qvalues = qvalues
                      )

    res <- res[ res$pvalue < pvalueCutoff, ]
    res <- res[ res$p.adjust < pvalueCutoff, ]
    idx <- order(res$pvalue, decreasing = FALSE)
    res <- res[idx, ]

    res$ID <- as.character(res$ID)
    row.names(res) <- res$ID

    if (verbose)
        print("done...")

    new("gseaResult",
        result     = res,
        setType    = setType,
        geneSets   = geneSets,
        geneList   = geneList,
        permScores = permScores,
        params     = params
        )
}


## GSEA algorithm (Subramanian et al. PNAS 2005)
## INPUTs to GSEA
## 1. Expression data set D with N genes and k samples.
## 2. Ranking procedure to produce Gene List L.
## Includes a correlation (or other ranking metric)
## and a phenotype or profile of interest C.
## 3. An exponent p to control the weight of the step.
## 4. Independently derived Gene Set S of N_H genes (e.g., a pathway).
## Enrichment Score ES.
## 2. Evaluate the fraction of genes in S ("hits") weighted
## by their correlation and the fraction of genes not in S ("miss")
## present up to a given position i in L.
gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
    ###################################################################
    ##    geneList                                                   ##
    ##                                                               ##
    ## 1. Rank order the N genes in D to form L = { g_1, ... , g_N}  ##
    ##    according to the correlation, r(g_j)=r_j,                  ##
    ##    of their expression profiles with C.                       ##
    ##                                                               ##
    ###################################################################

    ###################################################################
    ##    exponent                                                   ##
    ##                                                               ##
    ## An exponent p to control the weight of the step.              ##
    ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
    ##   the standard Kolmogorov-Smirnov statistic.                  ##
    ##   When p = 1, we are weighting the genes in S                 ##
    ##   by their correlation with C normalized                      ##
    ##   by the sum of the correlations over all of the genes in S.  ##
    ##                                                               ##
    ###################################################################

    ## genes defined in geneSet should appear in geneList.
    geneSet <- intersect(geneSet, names(geneList))

    N <- length(geneList)
    Nh <- length(geneSet)

    Phit <- Pmiss <- numeric(N)
    hits <- names(geneList) %in% geneSet ## logical

    Phit[hits] <- abs(geneList[hits])^exponent
    NR <- sum(Phit)
    Phit <- cumsum(Phit/NR)

    Pmiss[!hits] <-  1/(N-Nh)
    Pmiss <- cumsum(Pmiss)

    runningES <- Phit - Pmiss

    ## ES is the maximum deviation from zero of Phit-Pmiss
    max.ES <- max(runningES)
    min.ES <- min(runningES)
    if( abs(max.ES) > abs(min.ES) ) {
        ES <- max.ES
    } else {
        ES <- min.ES
    }

    if(fortify==TRUE) {
        df <- data.frame(x=seq_along(runningES),
                         runningScore=runningES,
                         position=as.integer(hits)
                         )
        return(df)
    }
    return (ES)
}

perm.geneList <- function(geneList) {
    perm.idx <- sample(seq_along(geneList), length(geneList), replace=FALSE)
    perm.geneList <- geneList
    names(perm.geneList) <- names(geneList)[perm.idx]
    return(perm.geneList)
}

perm.gseaEScore <- function(geneList, geneSet, nPerm, exponent=1) {
    res <- sapply(1:nPerm, function(i)
                  gseaScores(geneSet=geneSet,
                             geneList=perm.geneList(geneList),
                             exponent=exponent)
                  )
    return(res)
}
