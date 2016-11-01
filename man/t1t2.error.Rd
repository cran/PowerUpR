
\name{t1t2.error}
\alias{t1t2.error}
\title{Plots Type I and Type II Error Rates}

\description{
  \code{t1t2.error} plots Type I (\eqn{\alpha}) and Type II (\eqn{\beta}) error rates using central and
  noncentral \code{t} distributions for any objects returned from one of the \code{PowerUpR} functions.
}

\usage{
  t1t2.error(x)
}

\arguments{
  \item{x}{an object returned from one of the \code{PowerUpR} functions.}
}


\details{
\code{t1t2.error} plots Type I (\eqn{\alpha}) and Type II (\eqn{\beta}) error rates using central and
  noncentral \code{t} distributions for any objects returned from one of the \code{PowerUpR} functions.
  It is provided for didactic purpose.
}


\examples{
\dontrun{

     design1 <- mdes.bira2r1(rho2=.35, omega2=.10,
                            n=83, J=480)
     t1t2.error(design1)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}

\seealso{
  \code{\link{plot.pars}}
}

\keyword{plots}
