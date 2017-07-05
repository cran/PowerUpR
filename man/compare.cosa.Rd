
\name{compare.cosa}
\alias{compare.cosa}
\title{Compares COSA Solutions from Different Optimizers}

\description{
  \code{compare.cosa} compares COSA solutions from different optimizers for an object returned from one of the \code{cosa} functions.
}

\usage{
  compare.cosa(x)
}

\arguments{
  \item{x}{an object returned from one of the \code{cosa} functions.}
}


\details{
\code{compare.cosa} compares COSA solutions from different optimizers for an object returned from one of the \code{cosa} functions by passing parameters through \code{cosa} function.

}


\examples{
\dontrun{

     # object returned from cosa function
     design1 <- cosa.bira2r1(cn=1, cJ=10, cost=560,
                               constrain="cost", rho2=.20, omega2=.50)
     # compare COSA solution from different optimizers
     compare.cosa(design1)

  }
}

\keyword{object conversion}
