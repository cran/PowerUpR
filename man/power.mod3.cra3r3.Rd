\name{power.mod3.cra3r3}
\alias{power.mod3.cra3r3}
\title{Statistical Power Calculator for Moderator Effect in Three-Level Cluster Random Assignment Designs, Treatment and Moderator at Level-3}

\description{
  \code{power.mod3.cra3r3} calculates statistical power for moderator effect in designs with three-levels where level-3 units are randomly assigned to treatment and control groups and moderator is at level-3.
}

\usage{
  power.mod3.cra3r3(mdes=.25, alpha=.05, two.tail=TRUE,
                    rho2, rho3, P=.50, Q=NULL, g3=0, R12=0, R22=0, R32=0,
                    n, J, K, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{P}{proportion of level-3 units randomly assigned to treatment.}
  \item{Q}{proportion of level-1 units within moderator subgroup. \code{NULL} implies continuous moderator.}
  \item{g3}{number of covariates at level-3 excluding moderator.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{R32}{proportion of level-3 variance in the outcome explained by level-3 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{harmonic mean of level-2 units across level-3 units (or simple average).}
  \item{K}{level-3 sample size.}
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

     power.mod3.cra3r3(rho3=.1, rho2=.1, omega3=.05, Q=.5,
                    g3=1, R12=.3, R22=.4, R32=.5, n=20, J=4, K=60)

  }
}

\seealso{
  \code{\link{mdes.mod3.cra3r3}, \link{optimal.mod3.cra3r3}}
}

\keyword{power}
\keyword{moderator effect}
\keyword{cluster random assignment}
