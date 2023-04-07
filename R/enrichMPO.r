##' Enrichment analysis based on the DisGeNET (\url{http://www.disgenet.org/})
##'
##' given a vector of genes, this function will return the enrichment NCG
##' categories with FDR control
##'
##'
##' @inheritParams enrichNCG
##' @return A \code{enrichResult} instance
##' @export
##' @references Janet et al. (2015) DisGeNET: a discovery platform for the dynamical exploration of human diseases and their genes. \emph{Database} bav028
##' \url{http://database.oxfordjournals.org/content/2015/bav028.long}
##' @author Erqiang Hu
enrichMPO <- function(gene,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      universe,
                      minGSSize = 10,
                      maxGSSize = 500,
                      qvalueCutoff = 0.2,
                      readable = FALSE){

    enrichDisease(gene = gene,
                  pvalueCutoff = pvalueCutoff,
                  pAdjustMethod = pAdjustMethod,
                  universe = universe,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  qvalueCutoff = qvalueCutoff,
                  readable = readable,
                  ontology = "MPO")
    
}



#' Get MPO data
#' @param ont one of "DO" and "MPO".
#' @importFrom MPO.db MPOMGIDO
#' @importFrom MPO.db MPOANCESTOR
#' @importFrom MPO.db MPOPARENTS
#' @importFrom MPO.db MPOMPMGI
#' @importFrom MPO.db MPOTERM
#' @noRd 
get_MPO_data <- function(ont = "DO") {
    if (!exists(".DOSEenv")) .initial()
    .DOSEEnv <- get(".DOSEEnv", envir = .GlobalEnv)
    if (!exists(".MPO_DOSE_Env", envir=.DOSEEnv)) {
        assign(".MPO_DOSE_Env", new.env(), envir = .DOSEEnv)
        .MPO_DOSE_Env <- get(".MPO_DOSE_Env", envir = .DOSEEnv)
        if (ont == "DO") {
            ## get EXTID2PATHID gene_doid
            # (1) get eg.do
            # colnames(anno_old) <- c("eg", "doid")
            eg.do <- toTable(MPOMGIDO)
            colnames(eg.do) <- c("eg", "doid")
            # (2) DOSE:::rebuildAnnoData.internal(eg.do)
            DO2EG <- with(eg.do, split(as.character(eg), as.character(doid)))
            DOTERMs <- names(as.list(HDOANCESTOR))
            idx <- names(DO2EG) %in% DOTERMs
            DO2EG <- DO2EG[idx]
            DO2EG <- lapply(DO2EG, function(i) unique(i))
            EG2DO <- with(eg.do, split(as.character(doid), as.character(eg)))
            EG2DO <- lapply(EG2DO, function(i) unique(i[ i %in% DOTERMs ]))
            i <- unlist(lapply(EG2DO, function(i) length(i) != 0))
            EG2DO <- EG2DO[i]        
            EG2ALLDO <- lapply(EG2DO,
                               function(i) {
                                   ans <- unlist(mget(i, HDOANCESTOR))
                                   ans <- ans[ !is.na(ans) ]
                                   ans <- c(i, ans)
                                   ans <- unique(ans)
                                   return(ans)
                               })      
            len <- lapply(EG2ALLDO,length)
            EG2ALLDO.df <- data.frame(EG=rep(names(EG2ALLDO), times=len),
                                      DO=unlist(EG2ALLDO))
            DO <- NULL ## satisfy code tools
            DO2ALLEG <- with(EG2ALLDO.df, split(as.character(EG), as.character(DO)))
            DO2ALLEG <- lapply(DO2ALLEG, unique)

            PATH2NAME.df <- toTable(HDOTERM)
            # PATH2NAME.df <- PATH2NAME.df[, c("do_id", "Term")]
            PATH2NAME.df <- PATH2NAME.df[, c("doid", "term")]
            PATH2NAME.df <- unique(PATH2NAME.df)
            PATH2NAME <- PATH2NAME.df[,2]
            names(PATH2NAME) <- PATH2NAME.df[,1]
            
        } else {
            ## get EXTID2PATHID gene_doid
            # (1) get eg.do
            # colnames(anno_old) <- c("eg", "doid")
            eg.do <- toTable(MPOMPMGI)[, c(2,1)]
            colnames(eg.do) <- c("eg", "doid")
            # (2) DOSE:::rebuildAnnoData.internal(eg.do)
            DO2EG <- with(eg.do, split(as.character(eg), as.character(doid)))
            DOTERMs <- names(as.list(MPOANCESTOR))
            idx <- names(DO2EG) %in% DOTERMs
            DO2EG <- DO2EG[idx]
            DO2EG <- lapply(DO2EG, function(i) unique(i))
            EG2DO <- with(eg.do, split(as.character(doid), as.character(eg)))
            EG2DO <- lapply(EG2DO, function(i) unique(i[ i %in% DOTERMs ]))
            i <- unlist(lapply(EG2DO, function(i) length(i) != 0))
            EG2DO <- EG2DO[i]        
            EG2ALLDO <- lapply(EG2DO,
                               function(i) {
                                   ans <- unlist(mget(i, MPOANCESTOR))
                                   ans <- ans[ !is.na(ans) ]
                                   ans <- c(i, ans)
                                   ans <- unique(ans)
                                   return(ans)
                               })      
            len <- lapply(EG2ALLDO,length)
            EG2ALLDO.df <- data.frame(EG=rep(names(EG2ALLDO), times=len),
                                      DO=unlist(EG2ALLDO))
            DO <- NULL ## satisfy code tools
            DO2ALLEG <- with(EG2ALLDO.df, split(as.character(EG), as.character(DO)))
            DO2ALLEG <- lapply(DO2ALLEG, unique)
            PATH2NAME.df <- toTable(MPOTERM)
            colnames(PATH2NAME.df) <- c("doid", "term")
            PATH2NAME.df <- unique(PATH2NAME.df)
            PATH2NAME <- PATH2NAME.df[,2]
            names(PATH2NAME) <- PATH2NAME.df[,1]        
        }
        assign("EXTID2PATHID", EG2ALLDO, envir = .MPO_DOSE_Env)
        assign("PATHID2EXTID", DO2ALLEG, envir = .MPO_DOSE_Env)
        assign("PATHID2NAME", PATH2NAME, envir = .MPO_DOSE_Env)
    }
    get(".MPO_DOSE_Env", envir = .DOSEEnv)
}

