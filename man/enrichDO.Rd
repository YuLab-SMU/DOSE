\name{enrichDO}
\alias{enrichDO}
\title{DO Enrichment Analysis of a gene set.}
\usage{
  enrichDO(gene, ont = "DOLite", pvalueCutoff = 0.05,
    qvalueCutoff = 1, readable = F)
}
\arguments{
  \item{gene}{a vector of entrez gene id.}

  \item{ont}{one of DO or DOLite.}

  \item{pvalueCutoff}{Cutoff value of pvalue.}

  \item{qvalueCutoff}{Cutoff value of qvalue.}

  \item{readable}{whether mapping gene ID to gene Name}
}
\value{
  A \code{enrichResult} instance.
}
\description{
  Given a vector of genes, this function will return the
  enrichment DO categories with FDR control.
}
\examples{
set.seed(123)
	data(EG2DO)
	gene = sample(names(EG2DO), 30)
	yy = enrichDO(gene, pvalueCutoff=0.05)
	summary(yy)
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}
\seealso{
  \code{\link{enrichResult-class}}
}
\keyword{manip}

