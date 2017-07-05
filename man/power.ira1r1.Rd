
\name{power.ira1r1}
\alias{power.ira1r1}
\title{Statistical Power Calculator for Individual Random Assignment Designs, Completely Randomized Controlled Trials}

\description{
  \code{power.ira1r1} calculates statistical power for completely randomized controlled trials
    where individuals are randomly assigned to treatment or control groups.
}

\usage{
  power.ira1r1(mdes=.25, alpha=.05, two.tail=TRUE,
               P=.50, g1=0, R12=0,
               n, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{P}{proportion of units randomly assigned to treatment.}
  \item{g1}{number of covariates.}
  \item{R12}{proportion of variance in the outcome explained by covariates.}
  \item{n}{sample size.}
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

   power.ira1r1(n=55)

  }
}

\seealso{
  \code{\link{mdes.ira1r1}, \link{cosa.ira1r1}}
}

\keyword{power}
\keyword{main effect}
\keyword{individual random assignment}
