#' Shared parameters for DOSE functions
#'
#' @param gene a vector of entrez gene id
#' @param organism one of "hsa" and "mmu"
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
#' @name dose_params
NULL

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
#' @param organism organism
#' @param ont ont
#' @return DO Terms
#' @importMethodsFrom AnnotationDbi get
#' @importMethodsFrom AnnotationDbi exists
#' @export
#' @author Guangchuang Yu \url{https://yulab-smu.top}
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

