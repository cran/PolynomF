\name{as.function.polynom}
\alias{as.function.polynom}
\alias{as.function.polylist}

\title{Coerce polynom or polylist objects to function. }
\description{
  Since \code{polynom} objects are represented as functions,
  \code{as.function.polynom} simply removes the \code{class} attribute.
  The function \code{as.function.polylist} produces allows all
  polynomials on the list to be evaluated simultaneously at the same
  argument value.
}
\usage{
\method{as.function}{polynom}(x, variable = "x", \dots)
\method{as.function}{polylist}(x, \dots)
}

\arguments{
  \item{x}{A \code{polynom} or \code{polylist} object. }
  \item{variable}{A character string giving the name to be used for the
    formal argument of the resulting function, (for \code{as.polynom}
    only).}
  \item{\dots}{Not currently used. }
}
\details{
  Since \code{polynom} objects are already stored as functions,
  \code{as.function.polynom} is rarely needed and may be replaced by
  \code{unclass}.  An explicit call to \code{as.function.polynom},
  however, will generate a slightly faster version of the function as an
  unrooled loop, and does allow the user to specify a different name for
  the formal argument to be used.

  Since arithmetic on polynomials is provided, both of these
  functions may be evaluated at a \code{polynom} object as well as a
  numeric argument.
}
\value{
  An \code{R} function implementing the evaluation.
}
\references{None }
\author{Bill Venables}





\examples{
x <- polynom()
H <- polylist(polynom(1), x)
for(j in 2:5)
    H[[j+1]] <- x*H[[j]] - (j-1)*H[[j-1]]
H5 <- as.function(H)
H5(-4:4)
H5(x+1)
}


\keyword{symbolmath }

