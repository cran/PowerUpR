
\name{power.to.mdes}
\alias{power.to.mdes}
\title{power to MDES}

\description{
  \code{power.to.mdes} converts an object returned from \code{power} function into an object returned from \code{mdes} function.
}

\usage{
  power.to.mdes(x)
}

\arguments{
  \item{x}{an object returned from one of the \code{power} functions.}
}


\details{
  \code{power.to.mdes} converts an object returned from \code{power} function into an object returned from \code{mdes} function by passing parameters through \code{mdes} function.

}


\examples{
  \dontrun{

    # object returned from power function
    design1 <- power.bira2r1(mdes=.15, rho2=.35, omega2=.10, n=83, J=10)
    # convert the object into an object returned from mdes function
    design2 <- power.to.mdes(design1)

  }
}

\seealso{
  \code{\link{optimal.to.mdes}}
}

\keyword{object conversion}
