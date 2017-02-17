
\name{optimal.ira1r1}
\alias{optimal.ira1r1}
\title{Model 1.0: COSA Solver for Individual Random Assignment Designs, Completely Randomized Controlled Trials}

\description{
  \code{optimal.ira1r1} finds constrained optimal sample allocation (COSA) solutions for completely randomized controlled trials
    where individuals are randomly assigned to treatment and control groups.
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  optimal.ira1r1(cn, cost=NULL, n=NULL,
                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                 N0=c(10), ncase=10, gm=10,
                 constrain="cost", optimizer="auglag_cobyla",
                 P=.50, g1=0, R12=0)
}

\arguments{
  \item{cn}{marginal cost per unit.}
  \item{cost}{total cost or budget.}
  \item{n}{included for consistency, it should remain \code{NULL}.}
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{N0}{starting values for n.}
  \item{ncase}{number of cases to show in the output.}
  \item{gm}{grid multiplier to increase the range of sample size search.}
  \item{constrain}{one of the followings can be constrained at a specified cost or value: \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{optimizer}{algorithm to find optimal sample size given total cost, power, or MDES. Available algorithms: \code{"auglag_cobyla"}, \code{"auglag_lbfgs"}, \code{"auglag_mma"}, or \code{"auglag_slsqp"}.}
  \item{P}{proportion of units randomly assigned to treatment.}
  \item{g1}{number of covariates.}
  \item{R12}{proportion of variance in the outcome explained by covariates.}
}


\details{

An optimization is not necessary because the relationship between contraints and optimal sample is straight forward multiplication or division.
  Therefore use of this function is not recommended. Nonetheless, this function is provided for consistency and convenience.

Further definition of design parameters can be found in Dong & Maynard (2013).

}

\value{
  \item{fun}{function name.}
  \item{par}{list of parameters used in the function.}
  \item{nloptr}{list of \code{nloptr} log and output.}
  \item{round.optim}{solution after rounding. MDES is calculated at the specified power (default .80), and power is calculated at the specified MDES (default .25).}
  \item{integer.optim}{best integer solutions around \code{round.optim} solution. MDES is calculated at the specified power (default .80), and power is calculated at the specified MDES (default .25).}
}


\examples{
\dontrun{

     optimal.ira1r1(cn=1, cost=560,
                    constrain="cost")

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
  \code{\link{mdes.ira1r1}, \link{power.ira1r1}, \link{mrss.ira1r1}}
}
\keyword{optimal}
