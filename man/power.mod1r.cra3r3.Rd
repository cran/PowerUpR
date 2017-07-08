
\name{power.mod1r.cra3r3}
\alias{power.mod1r.cra3r3}
\title{Statistical Power Calculator for Randomly Varying Moderator Effect in Three-Level Cluster Random Assignment Designs, Treatment at Level-3, Moderator at Level-1}

\description{
  \code{power.mod1r.cra3r3} calculates statistical power for randomly varying moderator effect in designs with three-levels where level-3 units are randomly assigned to treatment and control groups and moderator is at level-1.
}

\usage{
  power.mod1r.cra3r3(mdes=.25, alpha=.05, two.tail=TRUE,
                      rho2, rho3, omega2, omega3, P=.50, Q=NULL,
                      g1=0, R12=0, RT22=0, RT32=0,
                      n, J, K, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{omega2}{moderator effect heterogeneity (unconditional) across level-2 units.}
  \item{omega3}{moderator effect heterogeneity (unconditional) across level-3 units.}
  \item{P}{proportion of level-3 units randomly assigned to treatment.}
  \item{Q}{proportion of level-1 units within moderator subgroup. \code{NULL} implies continuous moderator.}
  \item{g1}{number of covariates at level-1 excluding moderator.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{RT22}{proportion of level-2 variance in the moderator effect explained by level-2 covariates.}
  \item{RT32}{proportion of level-3 variance in the moderator effect explained by level-3 covariates.}
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

    power.mod1r.cra3r3(rho3=.05, rho2=.12, omega2=.08, omega3=.07,
                      P=.4, Q=.7, g1=1, R12=.20, RT22=0, RT32=0,
                      n=20, J=4, K=60)

  }
}

\seealso{
  \code{\link{mdes.mod1r.cra3r3}, \link{cosa.mod1r.cra3r3}}
}

\keyword{power}
\keyword{moderator effect}
\keyword{cluster random assignment}
