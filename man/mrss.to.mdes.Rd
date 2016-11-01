
\name{mrss.to.mdes}
\alias{mrss.to.mdes}
\title{MRSS to MDES}

\description{
  \code{mrss.to.mdes} converts an object returned from \code{mrss} function into an object returned from \code{mdes} function.
}

\usage{
  mrss.to.mdes(design)
}

\arguments{
  \item{design}{an object returned from one of the \code{mrss} functions.}
}


\details{
\code{mrss.to.mdes} converts an object returned from \code{mrss} function into an object returned from \code{mdes} function by passing parameters through \code{mdes} function.

}


\examples{
\dontrun{

     # object returned from mrss function
     design1 <- mrss.bira2r1(rho2=.35, omega2=.10, n=83)
     # convert the object into an object returned from mdes function
     design2 <- mrss.to.mdes(design1)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}


\seealso{
  \code{\link{mrss.to.power}, \link{optimal.to.mdes}, \link{optimal.to.power}}
}
\keyword{conversion}
