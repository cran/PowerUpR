\name{bcra3r2}
\alias{bcra3r2}
\alias{mdes.bcra3r2}
\alias{power.bcra3r2}
\alias{mrss.bcra3r2}
\title{Three-Level Blocked Cluster-level Random Assignment Design, Treatment at Level 2}

\description{
  Use \code{mdes.bcra3r2()} to calculate the minimum detectable effect size, \code{power.bcra3r2()} to calculate the statistical power, and \code{mrss.bcra3r2()} to calculate the minimum required sample size.
}

\usage{
mdes.bcra3r2(power=.80, alpha=.05, two.tailed=TRUE,
             rho2, rho3, omega3, p=.50, g3=0, r21=0, r22=0, r2t3=0,
             n, J, K)

power.bcra3r2(es=.25, alpha=.05, two.tailed=TRUE,
              rho2, rho3,  omega3, p=.50, g3=0, r21=0, r22=0, r2t3=0,
              n, J, K)

mrss.bcra3r2(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
             n, J, K0=10, tol=.10,
             rho2, rho3,  omega3, p=.50, g3=0, r21=0, r22=0, r2t3=0)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome between level 2 units (unconditional ICC2).}
  \item{rho3}{proportion of variance in the outcome between level 3 units (unconditional ICC3).}
  \item{omega3}{treatment effect heterogeneity as ratio of treatment effect variance among level 3 units to the residual variance at level 3.}
  \item{p}{average proportion of level 2 units randomly assigned to treatment within level 3 units.}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r2t3}{proportion of treatment effect variance among level 3 units explained by level 3 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{number of level 3 units.}
  \item{K0}{starting value for \code{K}.}
  \item{tol}{tolerance to end iterative process for finding \code{K}.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{df}{degrees of freedom.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mdes}{minimum detectable effect size.}
  \item{K}{number of level 3 units.}
}


\examples{
# cross-checks
mdes.bcra3r2(rho3=.13, rho2=.10, omega3=.4,
             n=10, J=6, K=24)
power.bcra3r2(es = .246, rho3=.13, rho2=.10, omega3=.4,
              n=10, J=6, K=24)
mrss.bcra3r2(es = .246, rho3=.13, rho2=.10, omega3=.4,
             n=10, J=6)
}

\seealso{
  \code{\link[cosa]{cosa.bcrd3r2}}
}
