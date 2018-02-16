\name{bcra4r2}
\alias{bcra4r2}
\alias{mdes.bcra4r2}
\alias{power.bcra4r2}
\alias{mrss.bcra4r2}
\title{Four-Level Blocked Cluster-level Random Assignment Design, Treatment at Level 2}

\description{
  \code{mdes.bcra4r2} calculates minimum detectable effect size, \code{power.bcra4r2} calculates statistical power, \code{mrss.bcra4r2} calculates minimum required sample size.
}

\usage{
mdes.bcra4r2(power=.80, alpha=.05, two.tailed=TRUE,
             rho2, rho3, rho4, omega3, omega4,
             p=.50, r21=0, r22=0, r2t3=0, r2t4=0, g4=0,
             n, J, K, L, ...)

power.bcra4r2(es=.25, alpha=.05, two.tailed=TRUE,
              rho2, rho3, rho4, omega3, omega4,
              p=.50, r21=0, r22=0, r2t3=0, r2t4=0, g4=0,
              n, J, K, L, ...)

mrss.bcra4r2(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
             n, J, K, L0=10, tol=.10,
             rho2, rho3, rho4, omega3, omega4,
             p=.50, r21=0, r22=0, r2t3=0, r2t4=0, g4=0, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level 3 units.}
  \item{rho4}{proportion of variance in the outcome explained by level 4 units.}
  \item{omega3}{treatment effect heterogeneity as ratio of treatment effect variance among level 3 units to the residual variance at level 3.}
  \item{omega4}{treatment effect heterogeneity as ratio of treatment effect variance among level 4 units to the residual variance at level 4.}
  \item{p}{average proportion of level 2 units randomly assigned to treatment within level 3 units.}
  \item{g4}{number of covariates at level 4.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r2t3}{proportion of treatment effect variance among level 3 units explained by level 3 covariates.}
  \item{r2t4}{proportion of treatment effect variance among level 4 units explained by level 4 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{harmonic mean of level 3 units across level 4 units (or simple average).}
  \item{L}{number of level 4 units.}
  \item{L0}{starting value for \code{L}.}
  \item{tol}{tolerance to end iterative process for finding \code{L}.}
  \item{...}{to handle depreciated or defunct arguments.}
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
mdes.bcra4r2(rho4=.05, rho3=.15, rho2=.15,
             omega4=.50, omega3=.50, n=10, J=4, K=4, L=20)

power.bcra4r2(es = .206, rho4=.05, rho3=.15, rho2=.15,
              omega4=.50, omega3=.50, n=10, J=4, K=4, L=20)

mrss.bcra4r2(es = .206, rho4=.05, rho3=.15, rho2=.15,
             omega4=.50, omega3=.50, n=10, J=4, K=4)
}

\seealso{
  \code{\link[cosa]{cosa.bcrd4r2}}
}
