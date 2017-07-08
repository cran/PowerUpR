\name{cosa.bira2f1}
\alias{cosa.bira2f1}
\title{COSA Solver for Main Effect in Two-Level Fixed Effects Blocked Individual Random Assignment Designs, Treatment at Level-1}

\description{
  \code{cosa.bira2f1} finds constrained optimal sample allocation (COSA) solutions for main effect in designs with two-levels
    where level-1 units are randomly assigned to treatment and control groups within level-2 units (fixed blocks).
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  cosa.bira2f1(cn=0, cJ=0, cost=NULL,
                n=NULL, J=NULL, P=NULL,
                nJ0=c(10,10), P0=.50,
                constrain="power", optimizer="auglag_slsqp",
                gm=2, ncase=10,
                power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                g1=0, R12=0)
}

\arguments{
  \item{cn}{marginal cost per level-1 unit.}
  \item{cJ}{marginal cost per level-2 unit.}
  \item{cost}{total cost or budget.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{level-2 sample size.}
  \item{P}{average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{nJ0}{vector with a length of two to specifiy starting values for level-1 and level-2 sample sizes.}
  \item{P0}{starting value for average proportion of level-1 units randomly assigned to treatment within level-2 units.}
  \item{constrain}{one of the followings can be constrained at a specified cost or value: \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{optimizer}{algorithm to find optimal sample sizes given total cost, power, or MDES. Available algorithms: \code{"auglag_slsqp"}, \code{"auglag_lbfgs"}, \code{"auglag_mma"}  or \code{"auglag_cobyla"}.}
  \item{gm}{grid multiplier to increase the range of sample sizes for approximate integer solutions.}
  \item{ncase}{number of approximate integer solutions to be exported.}
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{g1}{number of covariates at level-1.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
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

     cosa.bira2f1(cn=1, cJ=10, cost=560,
                    constrain="cost")

  }
}

\seealso{
  \code{\link{mdes.bira2f1}, \link{power.bira2f1}}
}

\keyword{cosa}
\keyword{main effect}
\keyword{individual random assignment}
\keyword{block design}
