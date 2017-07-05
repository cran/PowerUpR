
\name{power.bira3r1}
\alias{power.bira3r1}
\title{Statistical Power Calculator for Three-Level Random Effects Blocked Individual Random Assignment Design, Individuals Randomized within Blocks}

\description{
  \code{power.bira3r1} calculates statistical power for designs with three-levels
    where level-1 units are randomly assigned to treatment and control groups within level-2 units (random blocks).
}

\usage{
  power.bira3r1(mdes=.25, alpha=.05, two.tail=TRUE,
               rho2, rho3, omega2, omega3,
               P=.50, R12=0, RT22=0, RT32=0, g3=0,
               n, J, K, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{omega2}{treatment effect heterogeneity as ratio of treatment effect variance among level-2 units to the residual variance at level-2.}
  \item{omega3}{treatment effect heterogeneity as ratio of treatment effect variance among level-3 units to the residual variance at level-3.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{g3}{number of covariates at level-3.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{RT22}{proportion of treatment effect variance among level-2 units explained by level-2 covariates.}
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
  \item{power}{statistical power \eqn{(1-\beta)}.}
}

\examples{
\dontrun{

    power.bira3r1(rho3=.20, rho2=.15, omega3=.10, omega2=.10,
                 n=69, J=10, K=100)
  }
}

\seealso{
  \code{\link{mdes.bira3r1}, \link{cosa.bira3r1}}
}

\keyword{power}
\keyword{main effect}
\keyword{individual random assignment}
\keyword{block design}
