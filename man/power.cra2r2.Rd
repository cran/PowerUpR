
\name{power.cra2r2}
\alias{power.cra2r2}
\title{Statistical Power Calculator for Main Effect in Two-Level Cluster Random Assignment Design, Treatment at Level-2}

\description{
  \code{power.cra2r2} calculates statistical power for main effect in designs with two-levels
    where level-2 units are randomly assigned to treatment and control groups.
}

\usage{
  power.cra2r2(mdes=.25, alpha=.05, two.tail=TRUE,
                     rho2, g2=0, P=.50, R12=0, R22=0,
                     n, J, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{P}{proportion of level-2 units randomly assigned to treatment.}
  \item{g2}{number of covariates at level-2.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{level-2 sample size.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
}

\examples{
\dontrun{

   power.cra2r2(rho2=.20, n=4, J=20)

  }
}

\seealso{
  \code{\link{mdes.cra2r2}, \link{cosa.cra2r2}}
}

\keyword{power}
\keyword{main effect}
\keyword{cluster random assignment}
