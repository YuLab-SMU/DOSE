\name{sim-methods}
\docType{methods}
\alias{sim}
\alias{sim-methods}
\alias{sim,DOParams-method}
\alias{DOParams}

\title{Methods for calculating semantic similarity}

\description{
  Given a \code{DOParams} instance containing a list of DO terms or Entrez Gene IDs, the ID type, the method for calculating similarity, and the combine method for combining similarity socres, this method will compute semantic similarity among DO terms or genes.
}

\arguments{
  \item{params}{A \code{DOParams} instance }
}
  
\examples{
	data(DO2EG)
	set.seed(123)
	id <- list(a=sample(names(DO2EG), 5),b= sample(names(DO2EG), 6))
	## Setting Parameters...
	params <- new("DOParams", IDs=id, type="DOID", method="Wang", combine="rcmax.avg")
	## Calculating Semantic Similarities...
	sim(params)
}

\seealso{
  \code{\link{DOParams-class}}
}
\keyword{methods}
