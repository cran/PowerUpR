
\name{power.crd2r2}
\alias{power.crd2r2}
\title{Statistical Power Calculator for Two-Level Cluster Regression Discontinuity Design, Treatment Based on Discontinuity at Level-2}

\description{
   Statistical power calculator for two-level cluster regression discontinuity design where treatment assignment is based on discontinuity at level-2.
}

\usage{
  power.crd2r2(mdes=.25, alpha=.05, two.tail=TRUE,
                           rho2, g2=0, P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                           R12=0, R22=0,
						   n, J, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{RTZ}{correlation between treatment status and assignment variable}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean).}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean).}
  \item{dist.Z}{distribution of assignment variable; either \code{normal} or \code{uniform}.}
  \item{g2}{number of covariates at level-2.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{number of level-2 units.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}}
 }

\examples{
\dontrun{

      power.crd2r2(rho2=.20, n=10, J=100)

  }
}

\seealso{
  \code{\link{cosa.crd2r2}, \link{mdes.crd2r2}}
}

\keyword{power}
\keyword{main effect}
\keyword{cluster regression discontinuity}