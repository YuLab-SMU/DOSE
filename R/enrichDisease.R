#' @importFrom enrichit setReadable
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

    organism <- match.arg(organism, c("hsa", "mm"))

    annoData <- get_anno_data(ontology)
    
    res <- enrichit:::enricher_internal(gene = gene,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe = universe,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             qvalueCutoff = qvalueCutoff,
                             gson = annoData)

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
        if (organism == "hsa") {
            res <- setReadable(res, 'org.Hs.eg.db')
        } else {
            res <- setReadable(res, 'org.Mm.eg.db')
        }
    }
    return(res)
}


get_anno_data <- function(ontology) {
    if (ontology == "NCG") {
        annoData <- get_NCG_data()
    } else if (ontology == "DisGeNET") {
        annoData <- get_DGN_data()
    } else if (ontology == "snpDisGeNET") {
        annoData <- get_VDGN_data()
    } else if (ontology %in% c("HDO", "MPO", "HPO")) {
        annoData <- get_dose_data(ontology)
    } else {
        stop("ontology not supported yet...")
    }
    
    return(annoData)
}

get_dose_data <- function(ontology = "HPO") {
    .DOSEEnv <- get_dose_env()
    .obj <- sprintf(".%s_DOSE_GSON", ontology)
    if (exists(.obj, envir=.DOSEEnv)) {
        res <- get(.obj, envir = .DOSEEnv)
        return(res)
    }

    TERM2ALLEG <- get_ont2allgene(ontology) 
    
    # Convert list to data.frame for gsid2gene
    gsid2gene <- stack(TERM2ALLEG)
    colnames(gsid2gene) <- c("gene", "gsid")
    gsid2gene <- gsid2gene[, c("gsid", "gene")]

    termmap <- GOSemSim:::get_onto_data(
        ontology, 
        table="term", 
        output = "data.frame")

    gsid2name <- unique(termmap[, c(1, 2)])
    colnames(gsid2name) <- c("gsid", "name")

    gson_obj <- gson::gson(gsid2gene = gsid2gene, 
                           gsid2name = gsid2name,
                           species = "Homo sapiens",
                           gsname = ontology,
                           keytype = "ENTREZID",
                           version = "unknown",
                           accessed_date = as.character(Sys.Date()))
    
    assign(.obj, gson_obj, envir = .DOSEEnv)
    return(gson_obj)    
}

