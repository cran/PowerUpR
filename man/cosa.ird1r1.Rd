
\name{cosa.ird1r1}
\alias{cosa.ird1r1}
\title{COSA for Main Effect in Individual Regression Discontinuity Design}

\description{
  \code{cosa.ird1r1} calculates constrained optimal sample for main effect in individual regression discontinuity design.
}

\usage{
  cosa.ird1r1(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
                       n0=10, tol=.10,
					   P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                       g1=0, R12=0)
}

\arguments{
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{n0}{starting value for sample size.}
  \item{tol}{tolerance.}
  \item{P}{average proportion of units randomly assigned to treatment.}
  \item{RTZ}{correlation between treatment status and assignment variable.}
  \item{k1}{left truncation point (in standard deviation untis from full normal distribution mean).}
  \item{k2}{right truncation point (in standard deviation untis from full normal distribution mean).}
  \item{dist.Z}{distribution of assignment variable; either \code{normal} or \code{uniform}.}
  \item{g1}{number of covariates.}
  \item{R12}{proportion of variance in the outcome explained by covariates.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in the function.}
  \item{n}{sample size.}
}

\examples{
\dontrun{

     cosa.ird1r1(mdes=.466)

  }
}

\seealso{
  \code{\link{mdes.ird1r1}, \link{power.ird1r1}}
}

\keyword{cosa}
\keyword{main effect}
\keyword{individual regression discontinuity}
