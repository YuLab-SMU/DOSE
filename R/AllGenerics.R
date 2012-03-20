##' @exportMethod sim
setGeneric(
           name = "sim",
           def=function(params){standardGeneric("sim")}
           )

##' @exportMethod "setReadable<-"
setGeneric(
           name="setReadable<-",
           def=function(x, value) {standardGeneric("setReadable<-")}
           )

##' Mapping External ID to Ontology Term ID
##'
##' @S3method EXTID2TERMID DO
##' @export
##' @param gene gene ID vector
##' @param organism organism
EXTID2TERMID <- function(gene, organism) {
    UseMethod("EXTID2TERMID")
}

##' Mapping Ontology Term ID to External ID
##'
##' @S3method TERMID2EXTID DO
##' @export
##' @param term term ID vector
##' @param organism organism
TERMID2EXTID <- function(term, organism) {
    UseMethod("TERMID2EXTID")
}

##' Get all background External ID.
##'
##' @S3method ALLEXTID DO
##' @export
##' @param organism organism
ALLEXTID <- function(organism) {
    UseMethod("ALLEXTID")
}

##' Mapping Ontology Term ID to Name Symbol or Description
##'
##' @S3method TERM2NAME DO
##' @export
##' @param term term ID vector
TERM2NAME <- function(term) {
    UseMethod("TERM2NAME")
}
