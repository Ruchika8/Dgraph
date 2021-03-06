% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/games.R
\name{sample_gnm}
\alias{gnm}
\alias{sample_gnm}
\title{Generate random graphs according to the G(n,m) Erdos-Renyi model}
\usage{
sample_gnm(n, m, directed = FALSE, loops = FALSE)

gnm(...)
}
\arguments{
\item{n}{The number of vertices in the graph.}

\item{m}{The number of edges in the graph.}

\item{directed}{Logical, whether the graph will be directed, defaults to
FALSE.}

\item{loops}{Logical, whether to add loop edges, defaults to FALSE.}

\item{...}{Passed to \code{sample_app}.}
}
\value{
A graph object.
}
\description{
This model is very simple, every possible edge is created with the same
constant probability.
}
\details{
The graph has \sQuote{n} vertices and \sQuote{m} edges,
and the \sQuote{m} edges are chosen uniformly randomly from the set of all
possible edges. This set includes loop edges as well if the \code{loops}
parameter is TRUE.
}
\examples{

g <- sample_gnm(1000, 1000)
degree_distribution(g)
}
\author{
Gabor Csardi \email{csardi.gabor@gmail.com}
}
\references{
Erdos, P. and Renyi, A., On random graphs, \emph{Publicationes
Mathematicae} 6, 290--297 (1959).
}
\seealso{
\code{\link{sample_gnp}}, \code{\link{sample_pa}}
}
\keyword{graphs}

