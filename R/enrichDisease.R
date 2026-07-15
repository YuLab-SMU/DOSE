#' @importFrom enrichit setReadable
#' @importFrom enrichit ora_gson
enrichDisease <- function(gene,
                          organism = NULL,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe,
                          minGSSize = 10,
                          maxGSSize = 500,
                          qvalueCutoff = 0.2,
                          readable = FALSE,
                          ontology){

    info <- .resolve_ontology_organism(ontology, organism)
    ontology <- info$ontology
    organism <- info$organism
    .validate_entrez_ids(gene)

    annoData <- get_anno_data(ontology)
    
    if (missing(universe)) universe <- NULL
    if (!is.null(universe)) .validate_entrez_ids(universe, "universe")
    res <- ora_gson(gene = gene,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe = universe,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             qvalueCutoff = qvalueCutoff,
                             gson = annoData)

    if (is.null(res))
        return(res)
    res@organism <- info$species
    
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
    } else if (ontology %in% c("HDO", "MPO", "HPO")) {
        annoData <- get_dose_data(ontology)
    } else {
        stop("ontology not supported yet...")
    }
    
    return(annoData)
}

get_dose_data <- function(ontology = "HPO") {
    info <- .resolve_ontology_organism(ontology)
    ontology <- info$ontology
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
                           species = info$species,
                           gsname = ontology,
                           keytype = "ENTREZID",
                           version = "unknown",
                           accessed_date = as.character(Sys.Date()))
    
    assign(.obj, gson_obj, envir = .DOSEEnv)
    return(gson_obj)    
}

