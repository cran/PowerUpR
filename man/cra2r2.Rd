\name{cra2r2}
\alias{cra2r2}
\alias{mdes.cra2r2}
\alias{power.cra2r2}
\alias{mrss.cra2r2}
\title{Two-Level Cluster-level Random Assignment Design, Treatment at Level 2}

\description{
  \code{mdes.cra2r2} calculates minimum detectable effect size, \code{power.cra2r2} calculates statistical power, \code{mrss.cra2r2} calculates minimum required sample size.
}

\usage{
mdes.cra2r2(power=.80, alpha=.05, two.tailed=TRUE,
            rho2, p=.50, g2=0, r21=0, r22=0,
            n, J, ...)

power.cra2r2(es=.25, alpha=.05, two.tailed=TRUE,
             rho2, g2=0, p=.50, r21=0, r22=0,
             n, J, ...)

mrss.cra2r2(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            n, J0=10, tol=.10,
            rho2, g2=0, p=.50, r21=0, r22=0, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{p}{proportion of level 2 units randomly assigned to treatment.}
  \item{g2}{number of covariates at level 2.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{level 2 sample size.}
  \item{J0}{starting value for \code{J}.}
  \item{tol}{tolerance to end iterative process for finding \code{J}.}
  \item{...}{to handle depreciated or defunct arguments.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{df}{degrees of freedom.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mdes}{minimum detectable effect size.}
  \item{J}{number of level 2 units.}
}


\examples{
# cross-checks
mdes.cra2r2(rho2=.17, n=15, J=20)

power.cra2r2(es=.629, rho2=.17, n=15, J=20)

mrss.cra2r2(es=.629, rho2=.17, n=15)
}

\seealso{
  \code{\link[cosa]{cosa.crd2r2}}
}
