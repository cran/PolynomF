\name{c.polylist}
\alias{c.polylist}
\alias{rep.polylist}
\alias{unique.polylist}
\title{Utility methods for manipulating polylist objects }
\description{
  These three functions implement methods for well-known generic
  functions. 
}
\usage{
\method{c}{polylist}(\dots, recursive = FALSE)
\method{rep}{polylist}(x, times, \dots)
\method{unique}{polylist}(x, incomparables = FALSE, \dots)
}

\arguments{
  \item{\dots}{As for the generic function. }
  \item{x}{As for the generic function.}
  \item{recursive}{As for the generic function.}
  \item{times}{As for the generic function.}
  \item{incomparables}{As for the generic function.}
}
\details{
  Perform familiar operations, retaining \code{polylist} class for the result.
}
\value{
  An object of class \code{polylist}
}
\references{None }
\author{Kurt Hornik, slightly modified by Bill Venables }

\examples{
p <- poly.from.zeros(-3:4)
p5 <- rep(polylist(p), 5)
p5
p6 <- c(p, p5)
p6
unique(p6)
}


\keyword{symbolmath }

