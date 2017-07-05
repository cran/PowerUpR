
\name{power.ird1r1}
\alias{power.ird1r1}
\title{Statistical Power Calculator for Main Effect in Individual Regression Discontinuity Design}

\description{
  \code{power.ird1r1} calculates statistical power for main effect in individual regression discontinuity design where treatment assignment is based on discontinuity at the individual level.
}

\usage{
  power.ird1r1(mdes=.25, alpha=.05, two.tail=TRUE,
                       P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                       g1=0, R12=0,
					   n, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{P}{average proportion of units randomly assigned to treatment.}
  \item{RTZ}{correlation between treatment status and assignment variable.}
  \item{k1}{left truncation point (in standard deviation units from full normal distribution mean).}
  \item{k2}{right truncation point (in standard deviation units from full normal distribution mean).}
  \item{dist.Z}{distribution of assignment variable; either \code{normal} or \code{uniform}.}
  \item{g1}{number of covariates.}
  \item{R12}{proportion of variance in the outcome explained by covariates.}
  \item{n}{sample size.}
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

     power.ird1r1(n=100)

  }
}

\seealso{
  \code{\link{cosa.ird1r1}, \link{mdes.ird1r1}}
}

\keyword{power}
\keyword{main effect}
\keyword{individual regression discontinuity}
