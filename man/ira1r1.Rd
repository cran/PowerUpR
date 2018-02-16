\name{ira1r1}
\alias{ira1r1}
\alias{mdes.ira1r1}
\alias{power.ira1r1}
\alias{mrss.ira1r1}
\title{Individual-level Random Assignment Design}

\description{
  \code{mdes.ira1r1} calculates minimum detectable effect size, \code{power.ira1r1} calculates statistical power, \code{mrss.ira1r1} calculates minimum required sample size.
}

\usage{
mdes.ira1r1(power=.80, alpha=.05, two.tailed=TRUE,
            p=.50, g1=0, r21=0, n, ...)

power.ira1r1(es=.25, alpha=.05, two.tailed=TRUE,
             p=.50, g1=0, r21=0, n, ...)

mrss.ira1r1(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            n0=10, tol=.10,
            p=.50, g1=0, r21=0, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{p}{proportion of units randomly assigned to treatment.}
  \item{g1}{number of covariates.}
  \item{r21}{proportion of variance in the outcome explained by covariates.}
  \item{n}{sample size.}
  \item{n0}{starting value for \code{n}.}
  \item{tol}{tolerance to end iterative process for finding \code{n}.}
  \item{...}{to handle depreciated or defunct arguments.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{df}{degrees of freedom.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mdes}{minimum detectable effect size.}
  \item{n}{sample size.}
}


\examples{
# cross-checks
mdes.ira1r1(n=250)

power.ira1r1(es=.356, n=250)

mrss.ira1r1(es=.356)
}

\seealso{
  \code{\link[cosa]{power.ird1r1}}
}
