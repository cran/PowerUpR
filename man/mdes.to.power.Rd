
\name{mdes.to.power}
\alias{mdes.to.power}
\title{MDES to power}

\description{
  \code{mdes.to.power} converts an object returned from \code{mdes} function into an object returned from \code{power} function.
}

\usage{
  mdes.to.power(x)
}

\arguments{
  \item{x}{an object returned from one of the \code{mdes} functions.}
}


\details{
  \code{mdes.to.power} converts an object returned from \code{mdes} function into an object returned from \code{power} function by passing parameters through \code{power} function.

}


\examples{
  \dontrun{

    # object returned from power function
    design1 <- mdes.bira2r1(rho2=.35, omega2=.10, n=83, J=10)
    # convert the object into an object returned from mdes function
    design2 <- mdes.to.power(design1)

  }
}

\seealso{
  \code{\link{optimal.to.mdes}, \link{power.to.mdes}}
}

\keyword{object conversion}
