\name{GCD}
\alias{GCD}
\alias{GCD.polynom}
\alias{GCD.polylist}
\alias{LCM}
\alias{LCM.polynom}
\alias{LCM.polylist}
\title{GCD and LCM or two or more polynomials }
\description{
  Functions to find the Greatest Common Divisor (GCD) or Least Common
  Multipe (LCM) of two or more polynomials, specified either as individual
  arguments or as a \code{polylist} object.
}
\usage{
\method{GCD}{polynom}(\dots)
\method{GCD}{polylist}(\dots)
\method{LCM}{polynom}(\dots)
\method{LCM}{polylist}(\dots)
}

\arguments{
  \item{\dots}{ Either individual \code{polynom} arguments or a single
    \code{polylist} object with all polynomials. }
}
\details{
  Uses the classical GCD and LCM algorithms with polynomial arithmetic.  
}
\value{
  A single polynomial object giving the GCD or LCM respectively,
  normalised to have the leading coefficient unity (i.e. a monic
  polynomial).
}
\references{ None. }
\author{ Kurt Hornik, slightly modified by Bill Venables. }

\examples{
p1 <- poly.from.zeros(-3:2)
p2 <- poly.from.zeros(0:4)

pgcd <- GCD(p1, p2)

pl <- polylist(p1, p2)
plcm <- LCM(pl)

polylist(pgcd, plcm)

}


\keyword{symbolmath}

