\name{bcra4f3}
\alias{bcra4f3}
\alias{mdes.bcra4f3}
\alias{power.bcra4f3}
\alias{mrss.bcra4f3}
\title{Four-Level (Fixed Treatment Effect) Blocked Cluster-level Random Assignment Design, Treatment at Level 3}

\description{
  Use \code{mdes.bcra4f3()} to calculate minimum detectable effect size, \code{power.bcra4f3()} to calculate statistical power, and \code{mrss.bcra4f3()} to calculate minimum required sample size.
}


\usage{
mdes.bcra4f3(power=.80, alpha=.05, two.tailed=TRUE,
             rho2, rho3, p=.50, r21=0, r22=0, r23=0, g3=0,
             n, J, K, L)

power.bcra4f3(es=.25, alpha=.05, two.tailed=TRUE,
              rho2, rho3, p=.50, r21=0, r22=0, r23=0, g3=0,
              n, J, K, L)

mrss.bcra4f3(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
             n, J, K, L0=10, tol=.10,
             rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome between level 2 units (unconditional ICC2).}
  \item{rho3}{proportion of variance in the outcome between level 3 units (unconditional ICC3).}
  \item{p}{average proportion of level 3 units randomly assigned to treatment within level 4 units.}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r23}{proportion of level 3 variance in the outcome explained by level 3 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{harmonic mean of level 3 units across level 4 units (or simple average).}
  \item{L}{number of level 4 units.}
  \item{L0}{starting value for \code{L}.}
  \item{tol}{tolerance to end iterative process for finding \code{L}.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{df}{degrees of freedom.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mdes}{minimum detectable effect size.}
  \item{L}{number of level 4 units.}
}

\examples{
# cross-checks
mdes.bcra4f3(rho3=.15, rho2=.15,
             n=10, J=4, K=4, L=15)
power.bcra4f3(es=0.339, rho3=.15, rho2=.15,
              n=10, J=4, K=4, L=15)
mrss.bcra4f3(es=0.339, rho3=.15, rho2=.15,
             n=10, J=4, K=4)
}
