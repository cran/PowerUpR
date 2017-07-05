
\name{cosa.to.power}
\alias{cosa.to.power}
\title{COSA to Power}

\description{
  \code{cosa.to.power} converts an object returned from \code{cosa} function into an object returned from \code{power} function.
}

\usage{
  cosa.to.power(x)
}

\arguments{
  \item{x}{an object returned from one of the \code{cosa} functions.}
}


\details{
\code{cosa.to.power} converts an object returned from \code{cosa} function into an object returned from \code{power} function by passing parameters through \code{power} function.

}


\examples{
\dontrun{

     # object returned from cosa function
     design1 <- cosa.bira2r1(cn=1, cJ=10, cost=560,
                               constrain="cost", rho2=.20, omega2=.50)
     # convert the object into an object returned from power function
     design2 <- cosa.to.power(design1)

  }
}

\seealso{
  \code{\link{cosa.to.mdes}}
}

\keyword{object conversion}
