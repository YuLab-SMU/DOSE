\name{enrichDOResult-class}
\docType{class}
\alias{enrichDOResult-class}
\alias{show,enrichDOResult-method}
\alias{summary,enrichDOResult-method}

\title{Class "enrichDOResult"}
\description{
  This class represents the result of DO enrichment analysis.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("enrichDOResult", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{enrichDOResult}:}{containing a dataframe of DO enrichment analysis result.}
	\item{\code{pvalueCutoff}:}{ Cutoff value of pvalue. }
	\item{\code{Organism}:}{ Organism. }
	\item{\code{Gene}:}{ input gene vector. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object="enrichDOResult")}}
	\item{summary}{\code{signature(object="enrichDOResult")}}
  }
}

\examples{
	set.seed(123)
	data(EG2DO)
	gene = sample(names(EG2DO), 30)
	yy = enrichDO(gene, pvalueCutoff=0.05)
	## yy is an enrichDOResult instance, which stored the enrichment analysis result.
	yy
	summary(yy)
}

\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\link{enrichDO}}
}

\keyword{classes}
