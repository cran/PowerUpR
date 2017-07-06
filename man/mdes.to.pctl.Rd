
\name{mdes.to.pctl}
\alias{mdes.to.pctl}
\title{MDES to Percentiles}

\description{
  \code{mdes.to.pctl} converts MDES values returned from \code{mdes} function into percentile values.
}

\usage{
  mdes.to.pctl(x)
}

\arguments{
  \item{x}{an object returned from \code{mdes} function or a numeric value in standard deviation units.}
}


\details{
  \code{mdes.to.pctl} converts MDES values returned from \code{mdes} function into percentile values. Values should be interpreted as percentile increase with reference to 50th percentile.

}


\examples{
  \dontrun{

    # MDES values returned from mdes function
    design1 <- mdes.cra3r3(power=.80, rho2=.06, rho3=.18,
                       g3=1, R12=.55, R22=.50, R32=.45,
                       P=.40, n=10, J=2, K=83)
    # convert MDES values into percentiles
    percentiles <- mdes.to.pctl(design1)

  }
}

\seealso{
  \code{\link{optimal.to.mdes}}
}

\keyword{object conversion}
