enrichDisease <- function(gene,
                          organism = "hsa",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe,
                          minGSSize = 10,
                          maxGSSize = 500,
                          qvalueCutoff = 0.2,
                          readable = FALSE,
                          ontology){

    if (ontology == "NCG") {
        annoData <- get_NCG_data()
    } else if (ontology == "DisGeNET") {
        annoData <- get_DGN_data()
    } else if (ontology == "snpDisGeNET") {
        annoData <- get_VDGN_data()
    } else if (ontology == "DO" || ontology == "DOLite") {
        if (organism == "hsa") {
            annoData <- get_DO_data(ontology)
        } else {
            annoData <- get_MPO_data(ont = "DO")
        }
        
    } else if (ontology == "MPO") {
        annoData <- get_MPO_data(ont = "MPO")
    } else if (ontology == "HPO") {
        annoData <- get_HPO_data()
    } else {
        stop("ontology not supported yet...")
    }
    
    res <- enricher_internal(gene = gene,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe = universe,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             qvalueCutoff = qvalueCutoff,
                             USER_DATA = annoData)

    if (is.null(res))
        return(res)
    if (organism == "hsa") {
        res@organism <- "Homo sapiens"
    } else {
        res@organism <- "Mus musculus"
    }
    
    res@keytype <- "ENTREZID"
    res@ontology <- ontology

    if(readable) {
        res <- setReadable(res, 'org.Hs.eg.db')
    }
    return(res)
}

##' @importFrom yulab.utils check_pkg
##' @importFrom yulab.utils get_fun_from_pkg
##' @importFrom utils stack
get_ont_info <- function(ontology) {
    ## selected columns of `genemap`
    cols <- c(2, 1)
    if (ontology == "HDO" || ontology == "DO") {
        check_pkg("HDO.db")
        # not used in current version
    } else if (ontology == "HPO") {
        check_pkg("HPO.db")
        genemap <- get_fun_from_pkg("HPO.db", "HPOGENE")
        ancmap <- get_fun_from_pkg("HPO.db", "HPOANCESTOR")
        termmap <- get_fun_from_pkg("HPO.db", "HPOTERM")
    } else if (ontology == "MPO") {
        check_pkg("MPO.db")
        genemap <- get_fun_from_pkg("MPO.db", "MPOMPMGI")
        ancmap <- get_fun_from_pkg("MPO.db", "MPOANCESTOR")
        termmap <- get_fun_from_pkg("MPO.db", "MPOTERM")
    } else if (ontology == "MDO") {
        check_pkg("MPO.db")
        genemap <- get_fun_from_pkg("MPO.db", "MPGMGIDO")
        cols <- c(1, 2)

        check_pkg("HPO.db")
        ancmap <- get_fun_from_pkg("HPO.db", "HPOANCESTOR")
        termmap <- get_fun_from_pkg("HPO.db", "HPOTERM")
    }
    # toTable(genemap)[, cols]
    res <- list(genemap = genemap,
            cols = cols,
            ancmap = ancmap,
            termmap = termmap
        )
}

get_dose_data <- function(ontology = "HPO") {
    if (!exists(".DOSEenv")) .initial()
    .DOSEEnv <- get(".DOSEEnv", envir = .GlobalEnv)
    .env <- sprintf(".%s_DOSE_Env", ontology)
    if (exists(.env, envir=.DOSEEnv)) {
        res <- get(.env, envir = .DOSEEnv)
        return(res)
    }

    assign(.env, new.env(), envir = .DOSEEnv)
    ret_env <- get(.env, envir = .DOSEEnv)

    ont_info <- get_ont_info(ontology)
    eg2term <- toTable(ont_info$genemap)[, ont_info$cols]
        
    TERMS <- names(as.list(ont_info$ancmap))
    i <- eg2term[,2] %in% TERMS
    eg2term <- eg2term[i,]
    # TERM2EG <- split(eg2term[,1], eg2term[,2])
    EG2TERM <- split(eg2term[,2], eg2term[,1])

    EG2ALLTERM <- lapply(EG2TERM,
                       function(i) {
                           ans <- unlist(mget(i, ont_info$ancmap))
                           ans <- ans[ !is.na(ans) ]
                           ans <- c(i, ans)
                           ans <- unique(ans)
                           return(ans)
                       })

    EG2ALLTERM.df <- unique(stack(EG2ALLTERM)[, c(2,1)])
    TERM2ALLEG <- split(EG2ALLTERM.df[,1], EG2ALLTERM.df[,2])

    PATH2NAME.df <- unique(toTable(ont_info$termmap))
    PATH2NAME <- setNames(PATH2NAME.df[,2], PATH2NAME.df[,1])        

    assign("EXTID2PATHID", EG2ALLTERM, envir = ret_env)
    assign("PATHID2EXTID", TERM2ALLEG, envir = ret_env)
    assign("PATHID2NAME", PATH2NAME, envir = ret_env)

    return(ret_env)    
}

