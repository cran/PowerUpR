
\name{cosa.mod2r.cra3r3}
\alias{cosa.mod2r.cra3r3}
\title{COSA Solver for Randomly Varying Moderator Effect in Three-Level Cluster Random Assignment Designs, Treatment at Level-3 and Moderator at Level-2}

\description{
  \code{cosa.mod2r.cra3r3} finds constrained optimal sample allocation (COSA) solutions for randomly varying moderator effect in designs with three-levels
    where level-3 units are randomly assigned to treatment and control groups and moderator is at level-2.
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  cosa.mod2r.cra3r3(cn=0, cJ=0, cK=0, cost=NULL,
                 n=NULL, J=NULL, K=NULL, P=NULL, nJK0=c(10,10,10), P0=.50,
				 constrain="power", optimizer="auglag_slsqp",
				 gm=2, ncase=10, 
                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                 rho2, rho3, omega3, Q=NULL,
                 g2=0, R12=0, R22=0, RT32=0)
}

\arguments{
  \item{cn}{marginal cost per level-1 unit.}
  \item{cJ}{marginal cost per level-2 unit.}
  \item{cK}{marginal cost per level-3 unit.}
  \item{cost}{total cost or budget.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{harmonic mean of level-2 units across level-3 units (or simple average).}
  \item{K}{level-3 sample size.}
  \item{P}{proportion of level-3 units randomly assigned to treatment.}
  \item{nJK0}{vector with a length of two to specifiy starting values for level-1, level-2 and level-3 sample sizes.}
  \item{P0}{starting value for proportion of level-3 units randomly assigned to treatment.}
  \item{constrain}{one of the followings can be constrained at a specified cost or value: \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{optimizer}{algorithm to find optimal sample sizes given total cost, power, or MDES. Available algorithms: \code{"auglag_cobyla"}, \code{"auglag_lbfgs"}, \code{"auglag_mma"}, or \code{"auglag_slsqp"}.}
  \item{gm}{grid multiplier to increase the range of sample sizes for approximate integer solutions.}
  \item{ncase}{number of approximate integer solutions to be exported.}
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{omega3}{moderator effect heterogeneity (unconditional) across level-3 units.}
  \item{Q}{proportion of level-1 units within moderator subgroup. \code{NULL} implies continuous moderator.}
  \item{g2}{number of covariates at level-2 excluding moderator.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{RT32}{proportion of level-3 variance in the moderator effect explained by level-3 covariates.}
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

      cosa.mod2r.cra3r3(cn=1, cJ=10, cK=100, n=10, constrain="mdes",
                            rho2=.20, rho3=.10, omega3=.10)

  }
}

\seealso{
  \code{\link{mdes.mod2r.cra3r3}, \link{power.mod2r.cra3r3}}
}
\keyword{cosa}
\keyword{moderator effect}
\keyword{cluster random assignment}
