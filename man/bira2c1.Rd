\name{bira2c1}
\alias{bira2c1}
\alias{mdes.bira2c1}
\alias{power.bira2c1}
\alias{mrss.bira2c1}
\title{Two-Level Blocked (Constant Treatment Effect) Individual-level Random Assignment Design, Treatment at Level 1}

\description{
  \code{mdes.bira2c1} calculates minimum detectable effect size, \code{power.bira2c1} calculates statistical power, \code{mrss.bira2c1} calculates minimum required sample size.
}

\usage{
  mdes.bira2c1(power=.80, alpha=.05, two.tailed=TRUE,
               p=.50, g1=0, r21=0,
               n, J, ...)

  power.bira2c1(es=.25, alpha=.05, two.tailed=TRUE,
               p=.50, g1=0, r21=0,
               n, J, ...)

  mrss.bira2c1(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
             n, J0=10, tol=.10,
             p=.50, g1=0, r21=0, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{p}{average proportion of level 1 units randomly assigned to treatment within level 2 units.}
  \item{g1}{number of covariates at level 1.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
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
mdes.bira2c1(n=15, J=20)

power.bira2c1(es=.325, n=15, J=20)

mrss.bira2c1(es=.325, n=15)
}
