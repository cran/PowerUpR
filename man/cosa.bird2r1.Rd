
\name{cosa.bird2r1}
\alias{cosa.bird2r1}
\title{COSA Solver for Two-Level Blocked Individual Regression Discontinuity Design, Treatment Based on Discontinuity at Level-1 (Random Treatment Effect)}

\description{
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit while minimizing sampling variance fo the treatment effect,
    (ii) under power constraints given marginal costs per unit while minimizing the total cost,
    (ii) under MDES constraints given marginal costs per unit while minimizing the total cost, and
    (iv) under sample size constraints for one or more levels along with any of the i ii, or iii options.
}

\usage{
  cosa.bird2r1(cn=0, cJ=0, cost=NULL,
	               n=NULL, J=NULL, nJ0=c(10,10),
	               constrain="power", optimizer="auglag_slsqp",
	               gm=2, ncase=10,
	               power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
	               rho2, omega2,
	               P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
	               g2=0, R12=0, RT22=0)
}

\arguments{
  \item{cn}{marginal cost per level-1 unit.}
  \item{cJ}{marginal cost per level-2 unit.}
  \item{cost}{total cost or budget.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{number of level-2 units.}
  \item{nJ0}{vector with a length of two to specifiy starting values for level-1, and level-2 sample sizes.}
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{constrain}{one of the followings can be constrained at a specified cost or value: \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{optimizer}{algorithm to find optimal sample sizes given total cost, power, or MDES. Available algorithms: \code{"auglag_slsqp"}, \code{"auglag_lbfgs"}, \code{"auglag_mma"}  or \code{"auglag_cobyla"}.}
  \item{gm}{grid multiplier to increase the range of sample sizes for approximate integer solutions.}
  \item{ncase}{number of approximate integer solutions to be exported.}
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{omega2}{treatment effect heterogeneity as ratio of treatment effect variance among level-2 units to the residual variance at level-2.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{RTZ}{correlation between treatment status and assignment variable}
  \item{k1}{left truncation point (in standard deviation untis from full normal distribution mean).}
  \item{k2}{right truncation point (in standard deviation untis from full normal distribution mean).}
  \item{dist.Z}{distribution of assignment variable; either \code{normal} or \code{uniform}.}
  \item{g2}{number of covariates at level-2.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{RT22}{proportion of treatment effect variance among level-2 units explained by level-2 covariates.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in the function.}
  \item{nloptr}{list of \code{nloptr} log and output.}
  \item{exact.optim}{exact solution. MDES is calculated at the specified power (default .80), and power is calculated at the specified MDES (default .25).}
  \item{round.optim}{rounded solution. MDES is calculated at the specified power (default .80), and power is calculated at the specified MDES (default .25).}
  \item{integer.optim}{approximate integer solutions. MDES is calculated at the specified power (default .80), and power is calculated at the specified MDES (default .25).}
}


\examples{
\dontrun{

     cosa.bird2r1(mdes=0.0446, cn=1, cJ=10, cost=560,
                    constrain="cost", rho2=.20, omega2=.50)

  }
}

\seealso{
  \code{\link{mdes.bird2r1}, \link{power.bird2r1}}
}

\keyword{cosa}
\keyword{main effect}
\keyword{individual regression discontinuity}
\keyword{block design}
