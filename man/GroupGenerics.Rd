% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolynomF.R
\name{GroupGenerics}
\alias{GroupGenerics}
\alias{Summary.polynom}
\alias{Summary.polylist}
\alias{Math.polynom}
\alias{Math.polylist}
\title{Summary and Math methods for polynomials}
\usage{
\method{Summary}{polynom}(..., na.rm = FALSE)

\method{Summary}{polylist}(..., na.rm = FALSE)

\method{Math}{polynom}(x, ...)

\method{Math}{polylist}(x, ...)
}
\arguments{
\item{...}{Additional arguments}

\item{na.rm}{Logical: should missing values be removed?}

\item{x}{a \code{"polynom"} or \code{"polylist"} objects.}
}
\value{
The result of the group generic operation
}
\description{
These provide methods for the generic function \code{Summary}
and \code{Math} for polynomial and polylist objects.  For \code{Summary}
only \code{sum} and \code{prod} members are implemented
}
\examples{
lis <- as_polylist(lapply(-2:3, function(x) polynomial() - x))
prod(lis)
sum(lis)
solve(prod(lis))
solve(sum(lis))
}
