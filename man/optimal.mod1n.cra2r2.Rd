\name{optimal.mod1n.cra2r2}
\alias{optimal.mod1n.cra2r2}
\title{COSA Solver for Non-randomly Varying Moderator Effect in Two-Level Cluster Random Assignment Design, Treatment at level-2 and Moderator at Level-1}

\description{
  \code{optimal.mod1n.cra2r2} finds constrained optimal sample allocation (COSA) solutions for non-randomly varying moderator effect in designs with two-levels
    where level-2 units are randomly assigned to treatment and control groups and moderator is at level-1.
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  optimal.mod1n.cra2r2(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                  power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                  nJ0=c(10,10), ncase=10, gm=2,
                  constrain="power", optimizer="auglag_cobyla",
                  rho2, P=.50, Q=NULL, g1=0, R12=0)
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
  \item{P}{proportion of level-2 units randomly assigned to treatment.}
  \item{Q}{proportion of level-1 units within moderator subgroup. \code{NULL} implies continuous moderator.}
  \item{g1}{number of covariates at level-1 excluding moderator.}
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

     optimal.mod1n.cra2r2(cn=1, cJ=10, cost=560,
                          constrain="cost", rho2=.20)

  }
}

\seealso{
  \code{\link{mdes.mod1n.cra2r2}, \link{power.mod1n.cra2r2}}
}

\keyword{optimal}
\keyword{moderator effect}
\keyword{cluster random assignment}