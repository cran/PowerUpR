
\name{mrss.cra3r3}
\alias{mrss.cra3r3}
\title{Model 3.2: MRSS Calculator for 3-Level Cluster Random Assignment Designs, Treatment at Level 3}

\description{
  \code{mrss.cra3r3} calculates minimum required sample size for designs with 3-levels
    where level 3 units are randomly assigned to treatment and control groups.
}

\usage{
  mrss.cra3r3(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
              gm=2, ncase=10, constrain="power",
              n=NULL, J=NULL, K=NULL, K0=10, tol=.10,
              rho2, rho3,
              P=.50, g3=0, R12=0, R22=0, R32=0)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{power}{statistical power (1 - type II error).}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{gm}{grid multiplier to increase the range of sample size search for each level.}
  \item{ncase}{number of cases to show in the output.}
  \item{constrain}{parameter to contrain; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{level 3 sample size.}
  \item{K0}{starting value for estimating number of level 3 units.}
  \item{tol}{tolerance to stop the search algorithm.}
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level 3 units.}
  \item{P}{proportion of level 3 units randomly assigned to treatment.}
  \item{g3}{number of covariates at level 3.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{R22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{R32}{proportion of level 3 variance in the outcome explained by level 3 covariates.}
}


\details{

Level 3 sample size (\code{K}) is calculated
using an iterative procedure described in Dong & Maynard (2013) due to degrees of freedom dependency on \code{K}.
For other levels (\code{n}, and \code{J}) MRSS calculation is simply solving for the unknown.
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

     mrss.cra3r3(rho3=.06, rho2=.17,
                 n=15, J=3)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}

\references{


Dong, N., & Maynard, R. A. (2013). PowerUp!: A Tool for Calculating Minum Detectable Effect Sizes and Minimum Required Sample Sizes for Experimental and Quasi-Experimental Design Studies,\emph{Journal of Research on Educational Effectiveness, 6(1)}, 24-6.


}

\seealso{
  \code{\link{mdes.cra3r3}, \link{power.cra3r3}, \link{optimal.cra3r3}}
}
\keyword{mrss}
