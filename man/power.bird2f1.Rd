
\name{power.bird2f1}
\alias{power.bird2f1}
\title{Statistical Power Calculator for Two-Level Blocked Individual Regression Discontinuity Design, Treatment Based on Discontinuity at Level-1 (Fixed Treatment Effect)}

\description{
   Statistical power calculator for two-level blocked individual regression discontinuity design where treatment assignment is based on discontinuity at level-1 (fixed treatment effect).
}

\usage{
  power.bird2f1(mdes=.25, alpha=.05, two.tail=TRUE,
                           P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
						   g1=0, R12=0, 
						   n, J, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{RTZ}{correlation between treatment status and assignment variable}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean).}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean).}
  \item{dist.Z}{distribution of assignment variable; either \code{normal} or \code{uniform}.}
  \item{g1}{number of covariates at level-1.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
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

     power.bird2f1(n=10, J=100)

  }
}

\seealso{
  \code{\link{cosa.bird2f1}, \link{mdes.bird2f1}}
}

\keyword{power}
\keyword{main effect}
\keyword{individual regression discontinuity}
\keyword{block design}