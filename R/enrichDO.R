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
		
	GeneRatio <- mdply(data.frame(a=k, b=n), .yPaste)
	GeneRatio <- GeneRatio[,3]
	BgRatio <- mdply(data.frame(a=M, b=N), .yPaste)
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
		Organism = organism,
		Gene = gene
	)
}


setMethod("show", signature(object="enrichDOResult"),
	function (object){
		Organism = object@Organism
		GeneNum = length(object@Gene)
		pvalueCutoff=object@pvalueCutoff
		cat (GeneNum, Organism, "Genes to KEGG test for over-representation.", "\n", "p value <", pvalueCutoff, "\n")
	}
)

setMethod("summary", signature(object="enrichDOResult"),
	function(object) {
		return(object@enrichDOResult)
	}
)
