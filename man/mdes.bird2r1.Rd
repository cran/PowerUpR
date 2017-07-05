
\name{mdes.bird2r1}
\alias{mdes.bird2r1}
\title{MDES Calculator for Two-Level Blocked Individual Regression Discontinuity Design, Treatment Based on Discontinuity at Level-1 (Random Treatment Effect)}

\description{
   Minimum detectable effect size calculator for two-level blocked individual regression discontinuity design where treatment assignment is based on discontinuity at level-1 (random treatment effect).
}

\usage{
  mdes.bird2r1(power=.80, alpha=.05, two.tail=TRUE,
				rho2, omega2,
				P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
				g2=0, R12=0, RT22=0, n, J, ...)
}
				
\arguments{
  \item{power}{statistical power (1 - type II error).}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{omega2}{treatment effect heterogeneity as ratio of treatment effect variance among level-2 units to the residual variance at level-2.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{RTZ}{correlation between treatment status and assignment variable}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean).}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean).}
  \item{dist.Z}{distribution of assignment variable; either \code{normal} or \code{uniform}.}
  \item{g2}{number of covariates at level-2.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{RT22}{proportion of treatment effect variance among level-2 units explained by level-2 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{number of level-2 units.}
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

     mdes.bird2r1(rho2=.20, omega2=.50, n=10, J=100)

  }
}

\seealso{
  \code{\link{cosa.bird2r1}, \link{power.bird2r1}}
}

\keyword{mdes}
\keyword{main effect}
\keyword{individual regression discontinuity}
\keyword{block design}
