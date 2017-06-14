
\name{power.bcra3r2}
\alias{power.bcra3r2}
\title{Statistical Power Calculator for Three-Level Random Effects Blocked Cluster Random Assignment Designs, Treatment at Level-2}

\description{
  \code{power.bcra3r2} calculates statistical power for designs with three-levels
    where level-2 units are randomly assigned to treatment and control groups within level-3 units (random blocks).
}

\usage{
  power.bcra3r2(mdes=.25, alpha=.05, two.tail=TRUE,
       rho2, rho3, omega3,
       P=.50, g3=0, R12=0, R22=0, RT32=0,
       n, J, K, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{omega3}{treatment effect heterogeneity as ratio of treatment effect variance among level-3 units to the residual variance at level-3.}
  \item{P}{average proportion of level-2 units randomly assigned to treatment within level-3 units.}
  \item{g3}{number of covariates at level-3.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
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

    power.bcra3r2(rho3=.13, rho2=.10, omega3=.40,
                 n=10, J=6, K=24)

  }
}

\seealso{
  \code{\link{mdes.bcra3r2}, \link{optimal.bcra3r2}}
}

\keyword{power}
\keyword{main effect}
\keyword{cluster random assignment}
\keyword{block design}
