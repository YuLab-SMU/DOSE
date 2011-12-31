\name{infoContentMethod}
\alias{infoContentMethod}
\title{information content based methods}
\usage{
  infoContentMethod(ID1, ID2, ont = "DO", method,
    organism = "human")
}
\arguments{
  \item{ID1}{Ontology Term}

  \item{ID2}{Ontology Term}

  \item{ont}{Ontology}

  \item{method}{one of "Resnik", "Jiang", "Lin" and "Rel".}

  \item{organism}{one of supported species}
}
\value{
  semantic similarity score
}
\description{
  Information Content Based Methods for semantic similarity
  measuring
}
\details{
  implemented for methods proposed by Resnik, Jiang, Lin
  and Schlicker.
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}

