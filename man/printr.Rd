% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printr.R
\docType{package}
\name{printr}
\alias{printr}
\alias{printr-package}
\title{Create a printer callback function}
\format{An object of class \code{environment} of length 10.}
\usage{
printr
}
\arguments{
\item{fun}{The function to use as a printer callback function.}

\item{x}{An R object.}

\item{x}{The object to print, or a callback function. See
\code{\link{printer_callback}} for details.}

\item{max_lines}{Maximum number of lines to print, \emph{not}
including the header and the footer.}

\item{header}{The header, if a function, then it will be called,
otherwise printed using \code{cat}.}

\item{footer}{The footer, if a function, then it will be called,
otherwise printed using \code{cat}.}

\item{omitted_footer}{Footer that is only printed if anything
is omitted from the printout. If a function, then it will be called,
otherwise printed using \code{cat}.}

\item{...}{Extra arguments to pass to \code{print()}.}

\item{...}{Passed to the printing function.}

\item{.indent}{Character scalar, indent the printout with this.}

\item{.printer}{The printing function, defaults to \code{print}.}
}
\value{
\code{x}, invisibly.

The first element in \code{...}, invisibly.
}
\description{
A printer callback fucntion is a function can performs the actual
printing. It has a number of subcommands, that are called by
the \code{printer} package, in a form \preformatted{
    printer_callback("subcommand", argument1, argument2, ...)
} See the examples below.
}
\details{
The subcommands:

\describe{
  \item{\code{length}}{The length of the data to print, the number of
    items, in natural units. E.g. for a list of objects, it is the
    number of objects.}
  \item{\code{min_width}}{TODO}
  \item{\code{width}}{Width of one item, if \code{no} items will be
    printed. TODO}
  \item{\code{print}}{Argument: \code{no}. Do the actual printing,
    print \code{no} items.}
  \item{\code{done}}{TODO}
}
}
\keyword{datasets}

