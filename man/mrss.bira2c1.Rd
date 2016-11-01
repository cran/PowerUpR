
\name{mrss.bira2c1}
\alias{mrss.bira2c1}
\title{Model 2.1: MRSS Calculator for 2-Level Constant Effects Blocked Individual Random Assignment Designs, Treatment at Level 1}

\description{
  \code{mrss.bira2c1} calculates minimum required sample size (MRSS) for designs with 2-levels
    where level 1 units are randomly assigned to treatment and control groups within level 2 units (school intercepts only).
}

\usage{
  mrss.bira2c1(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
              gm=2, ncase=10, constrain="power",
              n=NULL, J=NULL, J0=10, n0=10, tol=.10,
              P=.50, g1=0, R12=0)
}
\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{power}{statistical power (1 - type II error).}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{gm}{grid multiplier to increase the range of sample size search for each level.}
  \item{ncase}{number of cases to show in the output.}
  \item{constrain}{parameter to constrain; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{level 2 sample size.}
  \item{n0}{starting value for estimating number of level 1 units.}
  \item{J0}{starting value for estimating number of level 2 units.}
  \item{tol}{tolerance to stop the search algorithm.}
  \item{P}{average proportion of level 1 units randomly assigned to treatment within level 2 units.}
  \item{g1}{number of covariates at level 1.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
}


\details{

Level 2 and level 1 sample sizes (\code{J} and \code{n}) are calculated
using an iterative procedure described in Dong & Maynard (2013) due to
model degrees of freedom dependency on \code{J} and \code{n}.
MRSS calculator returns values that are not integer. Rounding may produce
MDES and power values different from what was specified,
therefore an integer solution is approximated using brute force (See Value section).
Integer solution to MRSS for an omitted level assumes that specified sample sizes
for remaining levels may subject to some changes.

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

     mrss.bira2c1(n=83)

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
  \code{\link{mdes.bira2c1}, \link{power.bira2c1}, \link{optimal.bira2c1}}
}
\keyword{mrss}
