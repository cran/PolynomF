\name{as.character.polynom}
\alias{as.character.polynom}

\title{Coerce polynomial object to character}
\description{
  Coerces a \code{polynom} object to a printable character
  representation.
}
\usage{
\method{as.character}{polynom}(x, variable = "x", decreasing = FALSE, \dots)
}

\arguments{
  \item{x}{A \code{polynom} object}
  \item{variable}{Character string with the desired variable name. }
  \item{decreasing}{Should the powers be decreasing, or increasing as in
    power series form?}
  \item{\dots}{Not presently used }
}
\details{
  The character string may be parsed into an expression for the
  polynomial itself.
}
\value{
  A character string.
}
\references{ None }
\author{Bill Venables, with contributions from Kurt Hornik }

\examples{
p <- poly.from.zeros(-2:4)
as.character(p)
# "48*x - 28*x^2 - 56*x^3 + 35*x^4 + 7*x^5 - 7*x^6 + x^7"

}


\keyword{ symbolmath }

