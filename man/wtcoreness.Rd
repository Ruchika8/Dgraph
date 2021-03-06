% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structural.properties.R
\name{wtcoreness}
\alias{graph.wtcoreness}
\alias{wtcoreness}
\title{Weighted K-core decomposition of graphs}
\usage{
wtcoreness(graph, mode = c("all", "out", "in"), weights = NULL,
  type = c("first", "second"), alpha = NULL, beta = NULL, lambda = NULL,
  wtdeg = c("once", "rep"))
}
\arguments{
\item{graph}{The input graph, it can be directed or undirected}

\item{mode}{The type of the core in directed graphs. Character constant,
possible values: \code{in}: in-cores are computed, \code{out}: out-cores are
computed, \code{all}: the corresponding undirected graph is considered. This
argument is ignored for undirected graphs.}

\item{weights}{gives the weights of the edges in the graph. The order is determined 
by the edge ids. If not specified, the edge attribute weight will be used by default.}

\item{wtdeg}{The type of algorithm to be used. Character constant,
possible values: \code{once}: weighted degree is computed initially, i.e, once.
\code{rep}: weighted degree is computed initially as well as at each updation step.
If the graph is unweighted, the algorithm does not need this parameter.}
}
\value{
Numeric vector of integer numbers giving the weighted coreness of each
vertex.
}
\description{
The weighted k-core of a graph is the maximal subgraph in which each vertex has at least 
k as its weighted degree(combination of degree and weights of links incident to each node). 
The weighted coreness of a vertex is k if it belongs to the k-core but not to the (k+1)-core.
}
\details{
The cores of a graph form layers: the (k+1)-core is always a subgraph of the k-core.

This function calculates the weighted coreness for each vertex.
}
\examples{

g <- make_empty_graph(5)
g <- add_edges(g, c(1,2, 1,4, 1,5, 2,3, 2,5, 4,5))
x <- c(8,3,6,2,5,7)
wtcoreness(g,mode,x,wtdeg) 		

}
\author{
Chhaya Malik \email{chhaya.mcs.du.2015@gmail.com}

Mani Bansal \email{mani.mcs.du.2015@gmail.com}
}
\references{
Vladimir Batagelj, Matjaz Zaversnik: An O(m) Algorithm for Cores
Decomposition of Networks, 2002

Antonios Garas, Frank Schweitzer, and Shlomo Havlin: A k-shell decomposition method
for weighted networks. New Journal Of Physics, 2012.

Seidman S. B. (1983) Network structure and minimum degree, \emph{Social
Networks}, 5, 269--287.
}
\seealso{
\code{\link{wtcoreness}}
}
\keyword{graphs}

