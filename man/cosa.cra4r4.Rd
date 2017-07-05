
\name{cosa.cra4r4}
\alias{cosa.cra4r4}
\title{COSA Solver for Main Effect in Four-Level Cluster Random Assignment Designs, Treatment at Level-4}

\description{
  \code{cosa.cra4r4} finds constrained optimal sample allocation (COSA) solutions for main effect in designs with four-levels
    where level-4 units are randomly assigned to treatment and control groups.
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  cosa.cra4r4(cn=0, cJ=0, cK=0, cL=0, cost=NULL,
                n=NULL, J=NULL, K=NULL, L=NULL, P=NULL,
	              nJKL0=c(10,10,10,10), P0=.50,
                power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                rho4, rho3, rho2,
                g4=0, R42=0, R32=0, R22=0, R12=0)
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
  \item{P}{proportion of level-4 units randomly assigned to treatment.}
  \item{nJKL0}{vector with a length of four to specifiy starting values for level-1, level-2, level-3, and level-4 sample sizes.}
  \item{P0}{starting value for proportion of level-4 units randomly assigned to treatment.}
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{constrain}{one of the followings can be constrained at a specified cost or value: \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{optimizer}{algorithm to find optimal sample sizes given total cost, power, or MDES. Available algorithms: \code{"auglag_slsqp"}, \code{"auglag_lbfgs"}, \code{"auglag_mma"}  or \code{"auglag_cobyla"}.}
  \item{gm}{grid multiplier to increase the range of sample sizes for approximate integer solutions.}
  \item{ncase}{number of approximate integer solutions to be exported.}
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{rho4}{proportion of variance in the outcome explained by level-4 units.}
  \item{g4}{number of covariates at level-4.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{R32}{proportion of level-3 variance in the outcome explained by level-3 covariates.}
  \item{R42}{proportion of level-4 variance in the outcome explained by level-4 covariates.}
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

     cosa.cra4r4(cn=1, cJ=10, cK=100, cL=1000, cost=75600,
                    constrain="cost",
                    rho4=.05, rho3=.05, rho2=.10)

  }
}

\seealso{
  \code{\link{mdes.cra4r4}, \link{power.cra4r4}}
}
\keyword{cosa}
\keyword{main effect}
\keyword{cluster random assignment}
