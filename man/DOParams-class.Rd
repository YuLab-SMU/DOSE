\name{DOParams-class}
\docType{class}
\alias{DOParams-class}

\title{Class "DOParams"}
\description{
  A \code{DOParams} contains parameters for calculating DO semantic similarity among DO term or Gene list.
}

\section{Slots}{
  \describe{
    \item{\code{IDs}:}{containing a list of DO terms or Gene IDs.}
	\item{\code{type}:}{specify the type of IDs, one of "DOID", "GeneID".}
    \item{\code{ontology}:}{must set to "DO".}
	\item{\code{organism}:}{currently, only "human" supported.}
	\item{\code{method}:}{Method for calculating DO semantic similarity, one of "Resnik", "Jiang", "Lin", "Rel", "Wang".}
	\item{\code{combine}:}{Method for combining DO semantic similarity scores, one of "avg", "max", "rcmax", "rcmax.avg"}
	\item{\code{dropCodes}:}{dropCodes for mapping Gene to DO Terms.}	
  }
}

\section{Methods}{
  \describe{
    \item{sim}{
		\code{signature(object = "DOParams")}: calculates semantic similarities of \code{DOParams}.}
  }
}

\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
	\code{\link{sim}}
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
\keyword{classes}
