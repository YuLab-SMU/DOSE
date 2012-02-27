\name{enrich.internal}
\alias{enrich.internal}
\title{enrich.internal}
\usage{
  enrich.internal(gene, organism, pvalueCutoff,
    qvalueCutoff, ont, readable)
}
\arguments{
  \item{gene}{a vector of entrez gene id.}

  \item{organism}{supported organism.}

  \item{pvalueCutoff}{Cutoff value of pvalue.}

  \item{qvalueCutoff}{Cutoff value of qvalue.}

  \item{ont}{Ontology}

  \item{readable}{whether mapping gene ID to gene Name}
}
\value{
  A \code{enrichResult} instance.
}
\description{
  interal method for enrichment analysis
}
\details{
  using the hypergeometric model
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}
\keyword{manip}

