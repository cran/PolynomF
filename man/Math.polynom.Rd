\name{Math.polynom}
\alias{Math.polynom}
\alias{Math.polylist}
\alias{Ops.polynom}
\alias{Ops.polylist}
\alias{Summary.polynom}
\alias{Summary.polylist}

\title{S3 Group Generic Functions for Polynomial and Polylist Objects}
\description{
  The function \code{Ops.polynom} allows arithmetic operations on
  \code{polynom} objects: addition, subtraction, multiplication, division
  (with remainder), remainder, raising to a non-negaive integer power.  It
  also allows exact equality and exact inequality tests.

  \code{Math.polynom} allows \code{round}, \code{signif}, \code{floor},
  \code{ceiling} and \code{trunc} operations on a \code{polynom}
  coefficient vector, returning a \code{polynom} result.

  \code{Summary.polynom} and \code{Summary.polylist} allow \code{sum} and
  \code{prod} operations on \code{polynom} arguments, specified either as
  individual \code{polynom} objects or in a single \code{polylist} object.
}
\usage{
\method{Math}{polynom}(x, \dots)
\method{Math}{polylist}(x, \dots)
\method{Ops}{polynom}(e1, e2)
\method{Summary}{polynom}(\dots, na.rm = FALSE)
\method{Summary}{polylist}(\dots, na.rm = FALSE)
}

\arguments{
  \item{e1, e2}{Objects of class \code{polynom} or numeric arguments
    which may be coerced to class \code{polynom}.  At least one must be
    a \code{polynom} object.}
  \item{x}{An object of class \code{polymon} or \code{polylist}.}
  \item{\dots}{\code{polynom} objects or scalar numeric, as appropriate.
    In the case of \code{Summary.polylist}, may be a single
    \code{polylist} object.}
  \item{na.rm}{logical: should missing values be removed first?}
}
\details{
  None of these functions is called directly.  \code{Ops.polynom} is the
  workhorse of the entire package.  \code{Ops.polylist} allows
  vectorised arithmetic computations on \code{polylist}s.
}
\value{
  A \code{polynom} (or \code{polylist}) object giving the result of the operation.
}
\references{ None }
\author{Bill Venables, with some code provided by Kurt Hornik }





\examples{
x <- polynom(0:1)
p <- (x + 1)^2 - 3

round(p/3)
sum(p, p^2, p^3, p^4)

}


\keyword{symbolmath }

