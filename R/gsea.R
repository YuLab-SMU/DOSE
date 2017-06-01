##' @importFrom fgsea fgsea
GSEA_fgsea <- function(geneList,
                       exponent,
                       nPerm,
                       minGSSize,
                       maxGSSize,
                       pvalueCutoff,
                       pAdjustMethod,
                       verbose,
                       seed=FALSE,
                       USER_DATA) {

    if(verbose)
        message("preparing geneSet collections...")

    geneSets <- getGeneSet(USER_DATA)
    check_gene_id(geneList, geneSets)

    if(verbose)
        message("GSEA analysis...")

    tmp_res <- fgsea(pathways=geneSets,
                 stats=geneList,
                 nperm=nPerm,
                 minSize=minGSSize,
                 maxSize=maxGSSize,
                 gseaParam=exponent,
                 nproc = 0)

    p.adj <- p.adjust(tmp_res$pval, method=pAdjustMethod)
    qvalues <- calculate_qvalue(tmp_res$pval)

    Description <- TERM2NAME(tmp_res$pathway, USER_DATA)

    params <- list(pvalueCutoff = pvalueCutoff,
                   nPerm = nPerm,
                   pAdjustMethod = pAdjustMethod,
                   exponent = exponent,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize
                   )

    res <- data.frame(
        ID = as.character(tmp_res$pathway),
        Description = Description,
        setSize = tmp_res$size,
        enrichmentScore = tmp_res$ES,
        NES = tmp_res$NES,
        pvalue = tmp_res$pval,
        p.adjust = p.adj,
        qvalues = qvalues,
        stringsAsFactors = FALSE
    )

    res <- res[!is.na(res$pvalue),]
    res <- res[ res$pvalue <= pvalueCutoff, ]
    res <- res[ res$p.adjust <= pvalueCutoff, ]
    idx <- order(res$pvalue, decreasing = FALSE)
    res <- res[idx, ]

    if (nrow(res) == 0) {
        message("no term enriched under specific pvalueCutoff...")
        return(
            new("gseaResult",
                result     = res,
                geneSets   = geneSets,
                geneList   = geneList,
                params     = params,
                readable   = FALSE
                )
        )
    }

    row.names(res) <- res$ID
    observed_info <- lapply(geneSets[res$ID], function(gs)
        gseaScores(geneSet=gs,
                   geneList=geneList,
                   exponent=exponent)
        )

    if (verbose)
        message("leading edge analysis...")

    ledge <- leading_edge(observed_info)

    res$rank <- ledge$rank
    res$leading_edge <- ledge$leading_edge
    res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')

    if (verbose)
        message("done...")

    new("gseaResult",
        result     = res,
        geneSets   = geneSets,
        geneList   = geneList,
        params     = params,
        readable   = FALSE
        )
}

##' generic function for gene set enrichment analysis
##'
##'
##' @title GSEA_internal
##' @param geneList order ranked geneList
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param pvalueCutoff p value Cutoff
##' @param pAdjustMethod p value adjustment method
##' @param verbose print message or not
##' @param seed set seed inside the function to make result reproducible. FALSE by default.
##' @param USER_DATA annotation data
##' @param by one of 'fgsea' or 'DOSE'
##' @return gseaResult object
##' @author Yu Guangchuang
GSEA_internal <- function(geneList,
                 exponent,
                 nPerm,
                 minGSSize,
                 maxGSSize,
                 pvalueCutoff,
                 pAdjustMethod,
                 verbose,
                 seed=FALSE,
                 USER_DATA,
                 by="fgsea") {

    by <- match.arg(by, c("fgsea", "DOSE"))
    if (!is.sorted(geneList))
        stop("geneList should be a decreasing sorted vector...")
    if (by == 'fgsea') {
        .GSEA <- GSEA_fgsea
    } else {
        .GSEA <- GSEA_DOSE
    }

    res <- .GSEA(geneList          = geneList,
                 exponent          = exponent,
                 nPerm             = nPerm,
                 minGSSize         = minGSSize,
                 maxGSSize         = maxGSSize,
                 pvalueCutoff      = pvalueCutoff,
                 pAdjustMethod     = pAdjustMethod,
                 verbose           = verbose,
                 seed              = seed,
                 USER_DATA         = USER_DATA)
    res@organism <- "UNKNOWN"
    res@setType <- "UNKNOWN"
    res@keytype <- "UNKNOWN"
    return(res)
}

##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @importFrom stats p.adjust
##' @importFrom BiocParallel bplapply
##' @importFrom BiocParallel MulticoreParam
## @importFrom BiocParallel bpisup
## @importFrom BiocParallel bpstart
## @importFrom BiocParallel bpstop
##' @importFrom BiocParallel multicoreWorkers
GSEA_DOSE <- function(geneList,
                 exponent,
                 nPerm,
                 minGSSize,
                 maxGSSize,
                 pvalueCutoff,
                 pAdjustMethod,
                 verbose,
                 seed=FALSE,
                 USER_DATA) {

    if(verbose)
        message("preparing geneSet collections...")
    geneSets <- getGeneSet(USER_DATA)
    check_gene_id(geneList, geneSets)

    selected.gs <- geneSet_filter(geneSets, geneList, minGSSize, maxGSSize)

    if (is.null(selected.gs))
        return(NULL)


    if (verbose)
        message("calculating observed enrichment scores...")

    observed_info <- lapply(selected.gs, function(gs)
        gseaScores(geneSet=gs,
                   geneList=geneList,
                   exponent=exponent)
        )
    observedScore <- sapply(observed_info, function(x) x$ES)

    if (verbose) {
        message("calculating permutation scores...")
    }
    if (seed) {
        seeds <- sample.int(nPerm)
    }

    ## if (!bpisup()) {
    ##     bpstart(MulticoreParam(multicoreWorkers(), progressbar=verbose))
    ##     on.exit(bpstop())
    ## }

    permScores <- bplapply(1:nPerm, function(i) {
        if (seed)
            set.seed(seeds[i])
        perm.gseaEScore(geneList, selected.gs, exponent)
    })

    permScores <- do.call("cbind", permScores)

    rownames(permScores) <- names(selected.gs)

    pos.m <- apply(permScores, 1, function(x) mean(x[x >= 0]))
    neg.m <- apply(permScores, 1, function(x) abs(mean(x[x < 0])))


    normalized_ES <- function(ES, pos.m, neg.m) {
        s <- sign(ES)
        m <- numeric(length(ES))
        m[s==1] <- pos.m[s==1]
        m[s==-1] <- neg.m[s==-1]
        ES/m
    }

    NES <- normalized_ES(observedScore, pos.m, neg.m)

    permScores <- apply(permScores, 2, normalized_ES, pos.m=pos.m, neg.m=neg.m)

    if (verbose)
        message("calculating p values...")
    pvals <- sapply(seq_along(observedScore), function(i) {
        if( is.na(NES[i]) ) {
            NA
        } else if ( NES[i] >= 0 ) {
            (sum(permScores[i, ] >= NES[i]) +1) / (sum(permScores[i,] >= 0) +1)
        } else { # NES[i] < 0
            (sum(permScores[i, ] <= NES[i]) +1) / (sum(permScores[i,] < 0) +1)
        }

    })
    p.adj <- p.adjust(pvals, method=pAdjustMethod)
    qvalues <- calculate_qvalue(pvals)

    gs.name <- names(selected.gs)
    Description <- TERM2NAME(gs.name, USER_DATA)

    params <- list(pvalueCutoff = pvalueCutoff,
                   nPerm = nPerm,
                   pAdjustMethod = pAdjustMethod,
                   exponent = exponent,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize
                   )


    res <- data.frame(
        ID = as.character(gs.name),
        Description = Description,
        setSize = sapply(selected.gs, length),
        enrichmentScore = observedScore,
        NES = NES,
        pvalue = pvals,
        p.adjust = p.adj,
        qvalues = qvalues,
        stringsAsFactors = FALSE
    )

    res <- res[!is.na(res$pvalue),]
    res <- res[ res$pvalue <= pvalueCutoff, ]
    res <- res[ res$p.adjust <= pvalueCutoff, ]
    idx <- order(res$pvalue, decreasing = FALSE)
    res <- res[idx, ]

    if (nrow(res) == 0) {
        message("no term enriched under specific pvalueCutoff...")
        return(
            new("gseaResult",
                result     = res,
                geneSets   = geneSets,
                geneList   = geneList,
                params     = params,
                readable   = FALSE
                )
        )
    }

    row.names(res) <- res$ID
    observed_info <- observed_info[res$ID]

    if (verbose)
        message("leading edge analysis...")

    ledge <- leading_edge(observed_info)

    res$rank <- ledge$rank
    res$leading_edge <- ledge$leading_edge
    res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')


    if (verbose)
        message("done...")

    new("gseaResult",
        result     = res,
        geneSets   = geneSets,
        geneList   = geneList,
        permScores = permScores,
        params     = params,
        readable   = FALSE
        )
}


leading_edge <- function(observed_info) {
    core_enrichment <- lapply(observed_info, function(x) {
        runningES <- x$runningES
        runningES <- runningES[runningES$position == 1,]
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            leading_gene <- runningES$gene[1:i]
        } else {
            i <- which.min(runningES$runningScore)
            leading_gene <- runningES$gene[-c(1:i)]
        }
        return(leading_gene)
    })

    rank <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        ES <- x$ES
        if (ES >= 0) {
            rr <- which.max(runningES$runningScore)
        } else {
            i <- which.min(runningES$runningScore)
            rr <- nrow(runningES) - i + 1
        }
        return(rr)
    })

    tags <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        runningES <- runningES[runningES$position == 1,]
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            res <- i/nrow(runningES)
        } else {
            i <- which.min(runningES$runningScore)
            res <- (nrow(runningES) - i + 1)/nrow(runningES)
        }
        return(res)
    })

    ll <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            res <- i/nrow(runningES)
        } else {
            i <- which.min(runningES$runningScore)
            res <- (nrow(runningES) - i + 1)/nrow(runningES)
        }
        return(res)
    })

    N <- nrow(observed_info[[1]]$runningES)
    setSize <- sapply(observed_info, function(x) sum(x$runningES$position))
    signal <- tags * (1-ll) * (N / (N - setSize))

    tags <- paste0(round(tags * 100), "%")
    ll <- paste0(round(ll * 100), "%")
    signal <- paste0(round(signal * 100), "%")
    leading_edge <- paste0('tags=', tags, ", list=", ll, ", signal=", signal)

    res <- list(rank = rank,
                tags = tags,
                list = ll,
                signal = signal,
                leading_edge = leading_edge,
                core_enrichment = core_enrichment)
    return(res)
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
    ## geneSet <- intersect(geneSet, names(geneList))

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

    df <- data.frame(x=seq_along(runningES),
                     runningScore=runningES,
                     position=as.integer(hits)
                     )

    if(fortify==TRUE) {
        return(df)
    }

    df$gene = names(geneList)
    res <- list(ES=ES, runningES = df)
    return(res)
}

perm.geneList <- function(geneList) {
    ## perm.idx <- sample(seq_along(geneList), length(geneList), replace=FALSE)
    perm.idx <- sample.int(length(geneList))
    perm.geneList <- geneList
    names(perm.geneList) <- names(geneList)[perm.idx]
    return(perm.geneList)
}

perm.gseaEScore <- function(geneList, geneSets, exponent=1) {
    geneList <- perm.geneList(geneList)
    res <- sapply(1:length(geneSets), function(i)
                  gseaScores(geneSet=geneSets[[i]],
                             geneList=geneList,
                             exponent=exponent)$ES
                  )
    return(res)
}


geneSet_filter <- function(geneSets, geneList, minGSSize, maxGSSize) {
    geneSets <- sapply(geneSets, intersect, names(geneList))

    gs.idx <- get_geneSet_index(geneSets, minGSSize, maxGSSize)
    nGeneSet <- sum(gs.idx)

    if ( nGeneSet == 0 ) {
        msg <- paste0("No gene set have size between [", minGSSize, ", ", maxGSSize, "]...")
        message(msg)
        message("--> return NULL...")
        return(NULL)
    }
    geneSets[gs.idx]
}

