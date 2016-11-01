
\name{mrss.ira1r1}
\alias{mrss.ira1r1}
\title{Model 1.0: MRSS Calculator for Individual Random Assignment Designs, Completely Randomized Controlled Trials}

\description{
  \code{mrss.ira1r1} calculates minimum required sample size (MRSS) for completely randomized controlled trials
    where individuals are randomly assigned to treatment and control groups.
}

\usage{
  mrss.ira1r1(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
              gm=10, ncase=10, constrain="power", n=NULL, n0=10, tol=.10,
              P=.50, g1=0, R12=0)
}
\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{power}{statistical power (1 - type II error).}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{gm}{grid multiplier to increase the range of sample size search.}
  \item{ncase}{number of cases to show in the output.}
  \item{constrain}{parameter to constrain; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{n}{included for consistency, it should remain \code{NULL}.}
  \item{n0}{starting value for \code{n}}
  \item{tol}{tolerance to stop the search algorithm.}
  \item{P}{proportion of units randomly assigned to treatment.}
  \item{g1}{number of covariates.}
  \item{R12}{proportion of variance in the outcome explained by covariates.}
}


\details{

Sample size \code{n}) is calculated using an iterative procedure described in Dong & Maynard (2013)
due to model degrees of freedom dependency on \code{n}.

Further definition of design parameters can be found in Dong & Maynard (2013).

}

\value{
  \item{fun}{function name.}
  \item{par}{list of parameters used in MRSS calculation.}
  \item{round.mrss}{solution after rounding.}
  \item{integer.mrss}{best integer solutions around \code{round.mrss} solution.}
}


\examples{
\dontrun{

     mrss.ira1r1(n=83)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}

\references{


Dong & Maynard (2013). PowerUp!: A Tool for Calculating Minum Detectable Effect Sizes and Minimum Required Sample Sizes for Experimental and Quasi-Experimental Design Studies,\emph{Journal of Research on Educational Effectiveness, 6(1)}, 24-6.


}

\seealso{
  \code{\link{mdes.ira1r1}, \link{power.ira1r1}, \link{optimal.ira1r1}}
}
\keyword{mrss}
