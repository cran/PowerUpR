
\name{optimal.to.mdes}
\alias{optimal.to.mdes}
\title{Optimal to MDES}

\description{
  \code{optimal.to.mdes} converts an object returned from \code{optimal} function into an object returned from \code{mdes} function.
}

\usage{
  optimal.to.mdes(x)
}

\arguments{
  \item{x}{an object returned from one of the \code{optimal} functions.}
}


\details{
\code{optimal.to.mdes} converts an object returned from \code{optimal} function into an object returned from \code{mdes} function by passing parameters through \code{mdes} function.

}


\examples{
\dontrun{

     # object returned from optimal function
     design1 <- optimal.bira2r1(cn=1, cJ=10, cost=560,
                               constrain="cost", rho2=.20, omega2=.50)
     # convert the object into an object returned from mdes function
     design2 <- optimal.to.mdes(design1)

  }
}

\seealso{
  \code{\link{optimal.to.power}}
}

\keyword{object conversion}
