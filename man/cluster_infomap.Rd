% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/community.R
\name{cluster_infomap}
\alias{cluster_infomap}
\alias{infomap.community}
\title{Infomap community finding}
\usage{
cluster_infomap(graph, e.weights = NULL, v.weights = NULL, nb.trials = 10,
  modularity = TRUE)
}
\arguments{
\item{graph}{The input graph.}

\item{e.weights}{If not \code{NULL}, then a numeric vector of edge weights.
The length must match the number of edges in the graph.  By default the
\sQuote{\code{weight}} edge attribute is used as weights. If it is not
present, then all edges are considered to have the same weight.}

\item{v.weights}{If not \code{NULL}, then a numeric vector of vertex
weights. The length must match the number of vertices in the graph.  By
default the \sQuote{\code{weight}} vertex attribute is used as weights. If
it is not present, then all vertices are considered to have the same weight.}

\item{nb.trials}{The number of attempts to partition the network (can be any
integer value equal or larger than 1).}

\item{modularity}{Logical scalar, whether to calculate the modularity score
of the detected community structure.}
}
\value{
\code{cluster_infomap} returns a \code{\link{communities}} object,
please see the \code{\link{communities}} manual page for details.
}
\description{
Find community structure that minimizes the expected description length of a
random walker trajectory
}
\details{
Please see the details of this method in the references given below.
}
\examples{

## Zachary's karate club
g <- make_graph("Zachary")

imc <- cluster_infomap(g)
membership(imc)
communities(imc)

}
\author{
Martin Rosvall (\url{http://www.tp.umu.se/~rosvall/}) wrote the
original C++ code. This was ported to be more igraph-like by Emmanuel
Navarro (\url{http://enavarro.me/}).  The R interface and
some cosmetics was done by Gabor Csardi \email{csardi.gabor@gmail.com}.
}
\references{
The original paper: M. Rosvall and C. T. Bergstrom, Maps of
information flow reveal community structure in complex networks, \emph{PNAS}
105, 1118 (2008) \url{http://dx.doi.org/10.1073/pnas.0706851105},
\url{http://arxiv.org/abs/0707.0609}

A more detailed paper: M. Rosvall, D. Axelsson, and C. T. Bergstrom, The map
equation, \emph{Eur. Phys. J. Special Topics} 178, 13 (2009).
\url{http://dx.doi.org/10.1140/epjst/e2010-01179-1},
\url{http://arxiv.org/abs/0906.1405}.
}
\seealso{
Other community finding methods and \code{\link{communities}}.
}
\keyword{graphs}

