
\name{power.bcrd3r2}
\alias{power.bcrd3r2}
\title{Statistical Power Calculator for Three-Level Blocked Cluster Regression Discontinuity Design, Treatment Based on Discontinuity at Level-2 (Random Treatment Effect)}

\description{
     Statistical power calculator for three-level blocked cluster regression discontinuity design where treatment assignment is based on discontinuity at level-2 (random treatment effect).
}

\usage{
  power.bcrd3r2(mdes=.25, alpha=.05, two.tail=TRUE,
                           rho2, rho3, omega3,
                           RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                           P=.50, g3=0, R12=0, R22=0, RT32=0,
						   n, J, K, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{omega3}{treatment effect heterogeneity as ratio of treatment effect variance among level-3 units to the residual variance at level-3.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{RTZ}{correlation between treatment status and assignment variable}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean).}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean).}
  \item{dist.Z}{distribution of assignment variable; either \code{normal} or \code{uniform}.}
  \item{g3}{number of covariates at level-3.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{RT32}{proportion of treatment effect variance among level-3 units explained by level-3 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{harmonic mean of level-2 units across level-3 units (or simple average).}
  \item{K}{level-3 sample size.}
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

        power.bcrd3r2(rho3=.15, rho2=.20, omega3=.50, n=10, J=5, K=100)

  }
}
\seealso{
  \code{\link{cosa.bcrd3r2}, \link{mdes.bcrd3r2}}
}

\keyword{power}
\keyword{main effect}
\keyword{cluster regression discontinuity}
\keyword{block design}