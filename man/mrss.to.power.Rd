
\name{mrss.to.power}
\alias{mrss.to.power}
\title{MRSS to Power}

\description{
  \code{mrss.to.power} converts an object returned from \code{mrss} function into an object returned from \code{power} function.
}

\usage{
  mrss.to.power(design)
}

\arguments{
  \item{design}{an object returned from one of the \code{mrss} functions.}
}


\details{
\code{mrss.to.power} converts an object returned from \code{mrss} function into an object returned from \code{power} function by passing parameters through \code{power} function.

}


\examples{
\dontrun{

     # object returned from mrss function
     design1 <- mrss.bira2r1(rho2=.35, omega2=.10, n=83)
     # convert the object into an object returned from power function
     design2 <- mrss.to.power(design1)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}


\seealso{
  \code{\link{mrss.to.mdes}, \link{optimal.to.mdes}, \link{optimal.to.power}}
}
\keyword{conversion}
