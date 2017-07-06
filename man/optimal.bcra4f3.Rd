
\name{optimal.bcra4f3}
\alias{optimal.bcra4f3}
\title{COSA Solver for Main Effect in Four-Level Fixed Effects Blocked Cluster Random Assignment Designs, Treatment at Level-3}

\description{
  \code{optimal.bcra4f3} finds constrained optimal sample allocation (COSA) solutions for main effect in designs with four-levels
    where level-3 units are randomly assigned to treatment and control groups within level-4 units (fixed blocks).
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  optimal.bcra4f3(cn=0, cJ=0, cK=0, cL=0, cost=NULL,
                 n=NULL, J=NULL, K=NULL, L=NULL,
                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                 nJKL0=c(10,10,10,10), ncase=10, gm=2,
                 constrain="power", optimizer="auglag_cobyla",
                 rho3, rho2,
                 P=.50, g3=0, R32=0, R22=0, R12=0)
}

\arguments{
  \item{cn}{marginal cost per level-1 unit.}
  \item{cJ}{marginal cost per level-2 unit.}
  \item{cK}{marginal cost per level-3 unit.}
  \item{cL}{marginal cost per level-4 unit.}
  \item{cost}{total cost or budget.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{harmonic mean of level-2 units across level-3 units (or simple average).}
  \item{K}{harmonic mean of level-3 units across level-4 units (or simple average).}
  \item{L}{number of level-4 units.}
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{nJKL0}{vector with a length of four to specifiy starting values for level-1, level-2, level-3, and level-4 sample sizes.}
  \item{ncase}{number of cases to show in the output.}
  \item{gm}{grid multiplier to increase the range of sample size search for each level.}
  \item{constrain}{one of the followings can be constrained at a specified cost or value: \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{optimizer}{algorithm to find optimal sample sizes given total cost, power, or MDES. Available algorithms: \code{"auglag_cobyla"}, \code{"auglag_lbfgs"}, \code{"auglag_mma"}, or \code{"auglag_slsqp"}.}
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{P}{average proportion of level-3 units randomly assigned to treatment within level-4 units.}
  \item{g3}{number of covariates at level-3.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{R32}{proportion of level-3 variance in the outcome explained by level-3 covariates.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in the function.}
  \item{nloptr}{list of \code{nloptr} log and output.}
  \item{round.optim}{solution after rounding. MDES is calculated at the specified power (default .80), and power is calculated at the specified MDES (default .25).}
  \item{integer.optim}{best integer solutions around \code{round.optim} solution. MDES is calculated at the specified power (default .80), and power is calculated at the specified MDES (default .25).}
}

\examples{
\dontrun{

     optimal.bcra4f3(cn=1, cJ=10, cK=100, cL=1000, cost=75600,
                    constrain="cost",
                    rho3=.15, rho2=.20)

  }
}

\seealso{
  \code{\link{mdes.bcra4f3}, \link{power.bcra4f3}}
}

\keyword{optimal}
\keyword{main effect}
\keyword{cluster random assignment}
\keyword{block design}
