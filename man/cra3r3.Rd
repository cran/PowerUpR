\name{cra3r3}
\alias{cra3r3}
\alias{mdes.cra3r3}
\alias{power.cra3r3}
\alias{mrss.cra3r3}
\title{Three-Level Cluster-level Random Assignment Design, Treatment at Level 3}

\description{
  \code{mdes.cra3r3} calculates minimum detectable effect size, \code{power.cra3r3} calculates statistical power, \code{mrss.cra3r3} calculates minimum required sample size.
}

\usage{
mdes.cra3r3(power=.80, alpha=.05, two.tailed=TRUE,
            rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0,
            n, J, K, ...)

power.cra3r3(es=.25, alpha=.05, two.tailed=TRUE,
             rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0,
             n, J, K, ...)

mrss.cra3r3(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            n, J, K0=10, tol=.10,
            rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level 3 units.}
  \item{p}{proportion of level 3 units randomly assigned to treatment.}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r23}{proportion of level 3 variance in the outcome explained by level 3 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{level 3 sample size.}
  \item{K0}{starting value for \code{K}.}
  \item{tol}{tolerance to end iterative process for finding \code{K}.}
  \item{...}{to handle depreciated or defunct arguments.}
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
mdes.cra3r3(rho3=.06, rho2=.17, n=15, J=3, K=60)

power.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3, K=60)

mrss.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3)
}

\seealso{
  \code{\link[cosa]{cosa.crd3r3}}
}
