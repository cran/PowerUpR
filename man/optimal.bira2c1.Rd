\name{optimal.bira2c1}
\alias{optimal.bira2c1}
\title{COSA Solver for Main Effect in Two-Level Constant Effects Blocked Individual Random Assignment Designs, Treatment at Level-2}

\description{
  \code{optimal.bira2c1} finds constrained optimal sample allocation (COSA) solutions for main effct in designs with two-levels
    where level-1 units are randomly assigned to treatment and control groups within level-2 units (school intercepts only).
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  optimal.bira2c1(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                 nJ0=c(10,10), ncase=10, gm=2,
                 constrain="power", optimizer="auglag_cobyla",
                 P=.50, g1=0, R12=0)
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
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{g1}{number of covariates at level-1.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
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

     optimal.bira2c1(cn=1, cJ=10, cost=560,
                    constrain="cost")

  }
}

\seealso{
  \code{\link{mdes.bira2c1}, \link{power.bira2c1}}
}

\keyword{optimal}
\keyword{main effect}
\keyword{individual random assignment}
\keyword{block design}
