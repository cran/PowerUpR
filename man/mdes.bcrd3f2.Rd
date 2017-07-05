
\name{mdes.bcrd3f2}
\alias{mdes.bcrd3f2}
\title{MDES Calculator for Three-Level Blocked Cluster Regression Discontinuity Design,
Treatment Based on Discontinuity at Level-2 (Fixed Treatment Effect)}

\description{
  Minimum detectable effect size calculator for three-level blocked cluster regression discontinuity design where treatment assignment is based on discontinuity at level-2 (fixed treatment effect).
}

\usage{
  mdes.bcrd3f2(power=.80, alpha=.05, two.tail=TRUE,
                           rho2, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                           P=.50, g2=0, R12=0, R22=0,
						   n, J, K, ...)
}

\arguments{
  \item{power}{statistical power (1 - type II error).}
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
  \item{J}{harmonic mean of level-2 units across level-3 units (or simple average).}
  \item{K}{level-3 sample size.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in MDES calculation.}
  \item{ncp}{noncentrality parameter.}
  \item{mdes}{minimum detectable effect size and \eqn{(1-\alpha)}\% lower and upper confidence limits.}
 }


\examples{
	\dontrun{

      mdes.bcrd3f2(rho2=.20, n=10, J=5, K=100)

  }
}

\seealso{
  \code{\link{cosa.bcrd3f2}, \link{power.bcrd3f2}}
}

\keyword{mdes}
\keyword{main effect}
\keyword{cluster regression discontinuity}
\keyword{block design}
