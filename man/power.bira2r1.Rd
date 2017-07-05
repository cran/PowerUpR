\name{power.bira2r1}
\alias{power.bira2r1}
\title{Statistical Power Calculator for Two-Level Random Effects Blocked Individual Random Assignment Designs, Individuals Randomized within Blocks}

\description{
  \code{power.bira2r1} calculates statistical power for designs with two-levels
    where level-1 units are randomly assigned to treatment and control groups within level-2 units (random blocks).
}

\usage{
  power.bira2r1(mdes=.25, alpha=.05, two.tail=TRUE,
               rho2, omega2,
               g2=0, P=.50, R12=0, RT22=0,
               n, J, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{omega2}{treatment effect heterogeneity as ratio of treatment effect variance among level-2 units to the residual variance at level-2.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{g2}{number of covariates at level-2.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{RT22}{proportion of treatment effect variance among level-2 units explained by level-2 covariates.}
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

   power.bira2r1(rho2=.35, omega2=.10,
               n=83, J=480)

  }
}

\seealso{
  \code{\link{mdes.bira2r1}, \link{optimal.bira2r1}}
}

\keyword{power}
\keyword{main effect}
\keyword{individual random assignment}
\keyword{block design}