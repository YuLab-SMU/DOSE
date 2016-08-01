enrichDisease <- function(gene,
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
        annoData <- get_DO_data(ontology)
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
    
    res@organism <- "Homo sapiens"
    res@keytype <- "ENTREZID"
    res@ontology <- ontology

    if(readable) {
        res <- setReadable(res, 'org.Hs.eg.db')
    }
    return(res)
}
