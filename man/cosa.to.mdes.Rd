
\name{cosa.to.mdes}
\alias{cosa.to.mdes}
\title{COSA to MDES}

\description{
  \code{cosa.to.mdes} converts an object returned from \code{cosa} function into an object returned from \code{mdes} function.
}

\usage{
  cosa.to.mdes(x)
}

\arguments{
  \item{x}{an object returned from one of the \code{cosa} functions.}
}


\details{
\code{cosa.to.mdes} converts an object returned from \code{cosa} function into an object returned from \code{mdes} function by passing parameters through \code{mdes} function.

}


\examples{
\dontrun{

     # object returned from cosa function
     design1 <- cosa.bira2r1(cn=1, cJ=10, cost=560,
                               constrain="cost", rho2=.20, omega2=.50)
     # convert the object into an object returned from mdes function
     design2 <- cosa.to.mdes(design1)

  }
}

\seealso{
  \code{\link{cosa.to.power}}
}

\keyword{object conversion}
