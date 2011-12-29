##' Class "enrichDOResult"
##' This class represents the result of DO enrichment analysis.
##'
##'
##' @name enrichDOResult-class
##' @aliases enrichDOResult-class 
##'   show,enrichDOResult-method summary,enrichDOResult-method
##'
##' @docType class
##' @slot enrichDOResult DO enrichment result
##' @slot pvalueCutoff pvalueCutoff
##' @slot organism only "human" supported
##' @slot gene Gene IDs
##' @exportClass enrichDOResult
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{enrichDO}}
##' @keywords classes
setClass("enrichDOResult",
         representation=representation(
         enrichDOResult="data.frame",
         pvalueCutoff="numeric",
         organism = "character",
         gene = "character"
         )
         )


##' DO Enrichment Analysis of a gene set.
##'
##' Given a vector of genes, this function will return the enrichment DO
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param organism Currently, only "human" supported.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @return A \code{enrichDOResult} instance.
##' @importFrom plyr mdply
##' @importFrom qvalue qvalue
##' @importFrom methods new
##' @importFrom DO.db DOTERM
##' @importMethodsFrom DO.db Term
##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi mget
##' @importMethodsFrom AnnotationDbi exists
##' @importClassesFrom methods data.frame
##' @importClassesFrom methods list
##' @importClassesFrom DO.db DOTerms
##' @export
##' @seealso \code{\link{enrichDOResult-class}}
##' @author Guangchuang Yu \url{http://ygc.name}
##' @keywords manip
##' @examples
##'
##' 	set.seed(123)
##' 	data(EG2DO)
##' 	gene = sample(names(EG2DO), 30)
##' 	yy = enrichDO(gene, pvalueCutoff=0.05)
##' 	summary(yy)
##'
enrichDO <- function(gene, organism="human", pvalueCutoff=0.05) {
	if(!exists("DOSEEnv")) .initial()
	DO2EG <- get("DO2EG", envir=DOSEEnv)
	geneID.list = lapply(DO2EG, function(i) gene[gene %in% i])
	geneID <- sapply(geneID.list, function(i) paste(i, collapse="/"))
	k = sapply(geneID.list, length)
	M = sapply(DO2EG, length)
	doNum <- length(M)
	orgExtID <- unique(unlist(DO2EG))
	N <- rep(length(orgExtID), doNum)
	n <- rep(length(gene), doNum)

	args.df <- data.frame(numWdrawn=k-1, numW=M, numB=N-M, numDrawn=n)

	pvalues <- mdply(args.df, HyperG)
	pvalues <- pvalues[,5]

	GeneRatio <- mdply(data.frame(a=k, b=n), getRatio)
	GeneRatio <- GeneRatio[,3]
	BgRatio <- mdply(data.frame(a=M, b=N), getRatio)
	BgRatio <- BgRatio[,3]
	DOID <- names(DO2EG)
	term=mget(DOID, DOTERM, ifnotfound=NA)
	Description=sapply(term, Term)

	DOOver <- data.frame(DOID=DOID, Description=Description, GeneRatio=GeneRatio, BgRatio=BgRatio, pvalue=pvalues)

	#qvalue =  fdrtool(DOOver$pvalue, statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
	qobj = qvalue(DOOver$pvalue)
	qvalues <- qobj$qvalues
	DOOver <- data.frame(DOOver, qvalue=qvalues, geneID=geneID, Count=k)
	DOOver <- DOOver[order(pvalues),]
	DOOver <- DOOver[ DOOver$pvalue <= pvalueCutoff, ]
	DOOver <- DOOver[!is.na(DOOver$Description),]
	DOOver$Description <- as.character(DOOver$Description)

	new("enrichDOResult",
		enrichDOResult = DOOver,
		pvalueCutoff=pvalueCutoff,
		organism = organism,
		gene = gene
	)
}

##' show method for \code{enrichDOResult} instance
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @param object A \code{enrichDOResult} instance.
##' @return message
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("show", signature(object="enrichDOResult"),
	function (object){
		Organism = object@Organism
		GeneNum = length(object@Gene)
		pvalueCutoff=object@pvalueCutoff
		cat (GeneNum, Organism, "Genes to KEGG test for over-representation.", "\n", "p value <", pvalueCutoff, "\n")
	}
)

##' summary method for \code{enrichDOResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{enrichDOResult} instance.
##' @return A data frame
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("summary", signature(object="enrichDOResult"),
	function(object) {
		return(object@enrichDOResult)
	}
)
