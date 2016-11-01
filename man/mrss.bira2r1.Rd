
\name{mrss.bira2r1}
\alias{mrss.bira2r1}
\title{Model 2.3: MRSS Calculator for 2-Level Random Effects Blocked Individual Random Assignment Designs, Individuals Randomized within Blocks}

\description{
  \code{mrss.bira2r1} calculates minimum required sample size (MRSS) for designs with 2-levels
    where level 1 units are randomly assigned to treatment and control groups within level 2 units (random blocks).
}

\usage{
  mrss.bira2r1(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
              gm=2, ncase=10, constrain="power",
              n=NULL, J=NULL, J0=10, tol=.10,
              rho2, omega2,
              g2=0, P=.50, R12=0, RT22=0)
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
  \item{J0}{starting value for estimating number of level 2 units.}
  \item{tol}{tolerance to stop the search algorithm.}
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{omega2}{treatment effect heterogeneity as ratio of treatment effect variance among level 2 units to the residual variance at level 2.}
  \item{P}{average proportion of level 1 units randomly assigned to treatment within level 2 units.}
  \item{g2}{number of covariates at level 2.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{RT22}{proportion of treatment effect variance among level 2 units explained by level 2 covariates.}
}


\details{

Level 2 sample size (\code{J}) is calculated
using an iterative procedure described in Dong & Maynard (2013) due to
model degrees of freedom dependency on \code{J}.
For other level (\code{n}) MRSS calculation is simply solving for the unknown.
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

     mrss.bira2r1(rho2=.35, omega2=.10,
                 n=83)

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
  \code{\link{mdes.bira2r1}, \link{power.bira2r1}, \link{optimal.bira2r1}}
}
\keyword{mrss}
