
\name{optimal.to.power}
\alias{optimal.to.power}
\title{Optimal to Power}

\description{
  \code{optimal.to.power} converts an object returned from \code{optimal} function into an object returned from \code{power} function.
}

\usage{
  optimal.to.power(design)
}

\arguments{
  \item{design}{an object returned from one of the \code{optimal} functions.}
}


\details{
\code{optimal.to.power} converts an object returned from \code{optimal} function into an object returned from \code{power} function by passing parameters through \code{power} function.

}


\examples{
\dontrun{

     # object returned from optimal function
     design1 <- optimal.bira2r1(cn=1, cJ=10, cost=560,
                               constrain="cost", rho2=.20, omega2=.50)
     # convert the object into an object returned from power function
     design2 <- optimal.to.power(design1)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}


\seealso{
  \code{\link{optimal.to.mdes}, \link{mrss.to.mdes}, \link{mrss.to.power}}
}
\keyword{conversion}
