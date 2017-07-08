
\name{cosa.ira1r1}
\alias{cosa.ira1r1}
\title{Sample Size Calculator for Main Effect in Individual Random Assignment Designs, Completely Randomized Controlled Trials}

\description{
  \code{cosa.ira1r1} calculates minimum required sample size for main effect in completely randomized controlled trials where individuals are randomly assigned to treatment or control groups.
}

\usage{
  cosa.ira1r1(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
                       n0=10, tol=.10,
                       P=.50, g1=0, R12=0)
}

\arguments{
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{n0}{starting value for sample size.}
  \item{tol}{tolerance.}
  \item{P}{average proportion of units randomly assigned to treatment.}
  \item{g1}{number of covariates.}
  \item{R12}{proportion of variance in the outcome explained by covariates.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in the function.}
  \item{n}{sample size.}
}

\examples{
\dontrun{

     cosa.ira1r1(R12=.50)

  }
}

\seealso{
  \code{\link{mdes.ira1r1}, \link{power.ira1r1}}
}

\keyword{cosa}
\keyword{main effect}
\keyword{individual random assignment}
