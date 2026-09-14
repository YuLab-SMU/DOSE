#' Shared parameters for DOSE functions
#'
#' @param gene a vector of entrez gene id
#' @param organism species of the input Entrez gene IDs. Use "hsa" for human
#'   or "mm" for mouse. Common aliases such as "human", "mouse", and "mmu"
#'   are accepted. If omitted, the species is inferred from `ont`.
#' @param ont one of "HDO", "HPO" or "MPO"
#' @param pvalueCutoff pvalue cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes
#' @param minGSSize minimal size of genes annotated by ontology term for testing
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param qvalueCutoff qvalue cutoff
#' @param readable whether mapping gene ID to gene Name
#' @param geneList order ranked geneList
#' @param exponent weight of each step
#' @param nPerm permutation numbers
#' @param verbose print message or not
#' @param adaptive logical, use adaptive permutation or not (default: FALSE)
#' @param minPerm minimum number of permutations for adaptive mode (default: 1000)
#' @param maxPerm maximum number of permutations for adaptive mode (default: 10000)
#' @param method method of GSEA, one of "multilevel", "permute", "sample"
#' @param seed random seed for reproducibility, set to a number (or TRUE to use a
#'   fixed default seed) to make the result reproducible, or FALSE (default) to draw
#'   a random seed on each run, so results may vary between runs. The underlying
#'   permutation engine uses its own RNG seeded with this value; see
#'   \code{enrichit::gsea()} for details.
#' @name dose_params
NULL

.dose_ontology_info <- function(ontology) {
    if (length(ontology) != 1L || is.na(ontology) || !nzchar(ontology)) {
        stop("`ontology` must be one non-empty value.", call. = FALSE)
    }

    ontology <- toupper(ontology)
    if (ontology == "DO") ontology <- "HDO"

    organism <- c(HDO = "hsa", HPO = "hsa", MPO = "mm", NCG = "hsa")
    species <- c(
        HDO = "Homo sapiens",
        HPO = "Homo sapiens",
        MPO = "Mus musculus",
        NCG = "Homo sapiens"
    )

    if (!ontology %in% names(organism)) {
        stop(
            sprintf(
                "Unsupported ontology '%s'. Supported values are HDO, HPO, MPO, and NCG.",
                ontology
            ),
            call. = FALSE
        )
    }

    list(
        ontology = ontology,
        organism = unname(organism[[ontology]]),
        species = unname(species[[ontology]]),
        keytype = "ENTREZID"
    )
}

.normalize_organism <- function(organism) {
    if (length(organism) != 1L || is.na(organism) || !nzchar(organism)) {
        stop("`organism` must be one non-empty value.", call. = FALSE)
    }

    key <- tolower(organism)
    aliases <- c(
        hsa = "hsa",
        human = "hsa",
        `homo sapiens` = "hsa",
        mm = "mm",
        mmu = "mm",
        mouse = "mm",
        `mus musculus` = "mm"
    )

    if (!key %in% names(aliases)) {
        stop(
            sprintf(
                "Unsupported organism '%s'. Use 'hsa' for human or 'mm' for mouse.",
                organism
            ),
            call. = FALSE
        )
    }
    unname(aliases[[key]])
}

.resolve_ontology_organism <- function(ontology, organism = NULL) {
    info <- .dose_ontology_info(ontology)
    if (is.null(organism)) return(info)

    supplied <- .normalize_organism(organism)
    if (supplied != info$organism) {
        stop(
            sprintf(
                paste0(
                    "Ontology '%s' contains %s Entrez gene annotations and is ",
                    "incompatible with organism = '%s'. Use organism = '%s' or ",
                    "omit `organism` to infer it. Cross-species analysis requires ",
                    "an explicit ortholog conversion before enrichment."
                ),
                info$ontology,
                info$species,
                organism,
                info$organism
            ),
            call. = FALSE
        )
    }
    info
}

.validate_entrez_ids <- function(ids, arg = "gene") {
    ids <- as.character(ids)
    if (!length(ids)) {
        stop(sprintf("`%s` must contain at least one Entrez gene ID.", arg), call. = FALSE)
    }

    invalid <- is.na(ids) | !nzchar(ids) | !grepl("^[0-9]+$", ids)
    if (any(invalid)) {
        examples <- paste(utils::head(unique(ids[invalid]), 3L), collapse = ", ")
        stop(
            sprintf(
                "`%s` must contain Entrez gene IDs (`keytype = 'ENTREZID'`); invalid value(s): %s.",
                arg,
                examples
            ),
            call. = FALSE
        )
    }
    invisible(ids)
}

.validate_ranked_gene_list <- function(geneList) {
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("`geneList` must be a named numeric vector ranked by gene-level statistic.", call. = FALSE)
    }
    if (any(!is.finite(geneList))) {
        stop("`geneList` must contain only finite numeric values.", call. = FALSE)
    }
    .validate_entrez_ids(names(geneList), "names(geneList)")
    if (anyDuplicated(names(geneList))) {
        stop("`geneList` must have unique Entrez gene ID names.", call. = FALSE)
    }
    invisible(geneList)
}

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

#' @importFrom yulab.utils load_OrgDb
#' @importFrom GOSemSim get_organism
NULL



#' compute information content
#'
#'
#' @title compute information content
#' @param ont one of "DO", "HPO" and "MPO"
#' @return NULL
#' @importMethodsFrom AnnotationDbi toTable
#' @author Guangchuang Yu \url{https://yulab-smu.top}
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


#' provide gene ID, this function will convert to the corresponding DO Terms
#'
#'
#' @title convert Gene ID to DO Terms
#' @param gene entrez gene ID
#' @param organism species of the Entrez gene ID. If omitted, it is inferred
#'   from `ont`.
#' @param ont ont
#' @return DO Terms
#' @importMethodsFrom AnnotationDbi get
#' @importMethodsFrom AnnotationDbi exists
#' @export
#' @author Guangchuang Yu \url{https://yulab-smu.top}
gene2DO <- function(gene, organism = NULL, ont = "HDO") {
    info <- .resolve_ontology_organism(ont, organism)
    ont <- info$ontology
    .validate_entrez_ids(gene)
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

#' @importClassesFrom GOSemSim GOSemSimDATA
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

#' @importFrom memoise memoise
semdata2 <- memoise(semdata)


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

#' @importFrom ggplot2 facet_grid
#' @export
ggplot2::facet_grid


#' @importFrom GOSemSim get_organism
#' @export
GOSemSim::get_organism

#' @importFrom GOSemSim set_auto_update
#' @export
GOSemSim::set_auto_update

