setClass(Class="DOParams", 
	representation(
		IDs = "list",
		type = "character",
		method="character",
		combine="character",
		ontology="character", 
		dropCodes="character", 
		organism="character"),
	prototype=prototype (
		ontology="DO",
		dropCodes="NULL",
		organism="human")
	)
	

setClass("enrichDOResult",
	representation=representation(
		enrichDOResult="data.frame",
		pvalueCutoff="numeric",
		Organism = "character",
		Gene = "character"
	)
)	