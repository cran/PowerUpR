
\name{mdes.bira2f1}
\alias{mdes.bira2f1}
\title{MDES Calculator for Two-Level Fixed Effects Blocked Individual Random Assignment Designs, Treatment at Level-1}

\description{
  \code{mdes.bira2f1} calculates minimum detectable effect size (MDES) for designs with two-levels
    where level-1 units are randomly assigned to treatment and control groups within level-2 units (fixed blocks).
}

\usage{
  mdes.bira2f1(power=.80, alpha=.05, two.tail=TRUE,
              P=.50, g1=0, R12=0,
              n, J, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{g1}{number of covariates at level-1.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{level-2 sample size.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in MDES calculation.}
  \item{ncp}{noncentrality parameter.}
  \item{mdes}{minimum detectable effect size and \eqn{(1-\alpha)}\% lower and upper confidence limits.}
}

\examples{
\dontrun{

 mdes.bira2f1(n=55, J=3)

  }
}

\seealso{
  \code{\link{power.bira2f1}, \link{cosa.bira2f1}}
}

\keyword{mdes}
\keyword{main effect}
\keyword{individual random assignment}
\keyword{block design}
