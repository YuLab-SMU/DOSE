\name{wangMethod}
\alias{wangMethod}
\title{wangMethod}
\usage{
  wangMethod(ID1, ID2, ont = "DO", weight.isa = 0.8,
    weight.partof = 0.6, weight.do = 0.7)
}
\arguments{
  \item{ID1}{Ontology Term}

  \item{ID2}{Ontology Term}

  \item{ont}{Ontology}

  \item{weight.isa}{weight of isa relationship}

  \item{weight.partof}{weight of partof relationship}

  \item{weight.do}{weight of DO}
}
\value{
  semantic similarity score
}
\description{
  Method Wang for semantic similarity measuring
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}

