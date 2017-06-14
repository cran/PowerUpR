
\name{mdes.mod1r.cra2r2}
\alias{mdes.mod1r.cra2r2}
\title{MDES Calculator for Randomly Varying Moderator Effect in Two-Level Cluster Random Assignment Design, Treatment at Level-2, Moderator at Level-1}

\description{
  \code{mdes.mod1r.cra2r2} calculates minimum detectable effect size (MDES) for randomly varying moderator effect in designs with two-levels where level-2 units are randomly assigned to treatment and control groups and moderator is at level-1.
}

\usage{

  mdes.mod1r.cra2r2(power=.80, alpha=.05, two.tail=TRUE,
                     rho2, omega2, RT22=0, P=.50, Q=NULL, g1=0, R12=0,
                     n, J, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{omega2}{moderator effect heterogeneity (unconditional) across level-2 units.}
  \item{P}{proportion of level-2 units randomly assigned to treatment.}
  \item{Q}{proportion of level-1 units within moderator subgroup. \code{NULL} implies continuous moderator.}
  \item{g1}{number of covariates at level-1 excluding moderator.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{RT22}{proportion of level-2 variance in the moderator effect explained by level-2 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{level-2 sample size.}
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

 mdes.mod1r.cra2r2(rho2=.2, omega2=.2, RT22=.2, n=4, J=20)

  }
}

\seealso{
  \code{\link{power.mod1r.cra2r2}, \link{optimal.mod1r.cra2r2}}
}

\keyword{mdes}
\keyword{moderator effect}
\keyword{cluster random assignment}
