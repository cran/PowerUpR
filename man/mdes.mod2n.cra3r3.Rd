
\name{mdes.mod2n.cra3r3}
\alias{mdes.mod2n.cra3r3}
\title{MDES Calculator for Non-Randomly Varying Moderator Effect in Three-Level Cluster Random Assignment Design, Treatment at Level-3, Moderator at Level-2}

\description{
  \code{mdes.mod2n.cra3r3} calculates minimum detectable effect size (MDES) for non-randomly varying moderator effect in designs with three-levels where level-3 units are randomly assigned to treatment and control groups and moderator is at level-2.
}

\usage{
 mdes.mod2n.cra3r3(power=.80, alpha=.05, two.tail=TRUE,
                    rho2, rho3, P=.50, Q=NULL, g2=0, R12=0, R22=0,
                    n, J, K, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{P}{proportion of level-3 units randomly assigned to treatment.}
  \item{Q}{proportion of level-1 units within moderator subgroup. \code{NULL} implies continuous moderator.}
  \item{g2}{number of covariates at level-2 excluding moderator.}
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

  mdes.mod2n.cra3r3(rho3=.10, rho2=.10, Q=.5, g2=1,
                     R12=.30, R22=.40,  n=20, J=4, K=60)

  }
}

\seealso{
  \code{\link{power.mod2n.cra3r3}, \link{cosa.mod2n.cra3r3}}
}

\keyword{mdes}
\keyword{moderator effect}
\keyword{cluster random assignment}
