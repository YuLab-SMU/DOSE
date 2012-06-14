\name{geneSim}
\alias{geneSim}
\title{geneSim}
\usage{
  geneSim(geneID1, geneID2, method = "Wang",
    organism = "human", combine = "BMA")
}
\arguments{
  \item{geneID1}{entrez gene vector}

  \item{geneID2}{entrez gene vector}

  \item{method}{one of "Wang", "Resnik", "Rel", "Jiang",
  and "Lin".}

  \item{organism}{only "human" supported}

  \item{combine}{One of "max", "average", "rcmax", "BMA"
  methods, for combining semantic similarity scores of
  multiple DO terms associated with gene/protein.}
}
\value{
  score matrix
}
\description{
  measuring similarities bewteen two gene vectors.
}
\details{
  provide two entrez gene vectors, this function will
  calculate their similarity.
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}

