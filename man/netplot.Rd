\name{netplot}
\alias{netplot}
\title{netplot}
\usage{
  netplot(g, vertex.label.font = 2,
    vertex.label.color = "#666666", vertex.label.cex = 1.5,
    layout = layout.fruchterman.reingold,
    foldChange = NULL, fixed = TRUE, col.bin = 10)
}
\arguments{
  \item{g}{igraph object}

  \item{vertex.label.font}{font size}

  \item{vertex.label.color}{font text color}

  \item{layout}{layout}

  \item{foldChange}{fold change}

  \item{fixed}{logical}

  \item{col.bin}{number of legend color bin}
}
\value{
  plot
}
\description{
  plot network
}
\details{
  plot network of igraph object
}
\author{
  Yu Guangchuang
}

