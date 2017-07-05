
\name{optimal.bira2r1}
\alias{optimal.bira2r1}
\title{COSA Solver for Main Effect in Two-Level Random Effects Blocked Individual Random Assignment Designs, Individuals Randomized within Blocks}

\description{
  \code{optimal.bira2r1} finds constrained optimal sample allocation (COSA) solutions for main effect in designs with two-levels
    where level-1 units are randomly assigned to treatment and control groups within level-2 units (random blocks).
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  optimal.bira2r1(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                 nJ0=c(10,10), ncase=10, gm=2,
                 constrain="power", optimizer="auglag_cobyla",
                 rho2, omega2,
                 P=.50, g2=0, R12=0, RT22=0)
}

\arguments{
  \item{cn}{marginal cost per level-1 unit.}
  \item{cJ}{marginal cost per level-2 unit.}
  \item{cost}{total cost or budget.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{level-2 sample size.}
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{nJ0}{vector with a length of two to specifiy starting values for level-1 and level-2 sample sizes.}
  \item{ncase}{number of cases to show in the output.}
  \item{gm}{grid multiplier to increase the range of sample size search for each level.}
  \item{constrain}{one of the followings can be constrained at a specified cost or value: \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{optimizer}{algorithm to find optimal sample sizes given total cost, power, or MDES. Available algorithms: \code{"auglag_cobyla"}, \code{"auglag_lbfgs"}, \code{"auglag_mma"}, or \code{"auglag_slsqp"}.}
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{omega2}{treatment effect heterogeneity as ratio of treatment effect variance among level-2 units to the residual variance at level-2.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{g2}{number of covariates at level-2.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{RT22}{proportion of treatment effect variance among level-2 units explained by level-2 covariates.}
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

     optimal.bira2r1(cn=1, cJ=10, cost=560,
                    constrain="cost", rho2=.20, omega2=.50)

  }
}

\seealso{
  \code{\link{mdes.bira2r1}, \link{power.bira2r1}}
}

\keyword{optimal}
\keyword{main effect}
\keyword{individual random assignment}
\keyword{block design}
