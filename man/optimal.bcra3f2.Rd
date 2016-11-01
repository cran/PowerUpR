
\name{optimal.bcra3f2}
\alias{optimal.bcra3f2}
\title{Model 4.1: COSA Solver for 3-Level Fixed Effects Blocked Cluster Random Assignment Designs, Treatment at Level 2}

\description{
  \code{optimal.bcra3f2} finds constrained optimal sample allocation (COSA) solutions for designs with 3-levels
    where level 2 units are randomly assigned to treatment and control groups within level 3 units (fixed blocks).
    COSA can be found in the following forms,
    (i) under budgetary constraints given marginal costs per unit,
    (ii) under power contraints given marginal costs per unit,
    (iii) under MDES contraints given marginal costs per unit, and
    (iv) under sample size contraints for one or more levels along with any of the i,ii, or iii options.
}

\usage{
  optimal.bcra3f2(cn, cJ, cK, cost=NULL, n=NULL, J=NULL, K=NULL,
                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                 nJK0=c(10,10,10), ncase=10, gm=2,
                 constrain="cost", optimizer="auglag_cobyla",
                 rho2,
                 P=.50, g2=0,  R12=0, R22=0)
}

\arguments{
  \item{cn}{marginal cost per level 1 unit.}
  \item{cJ}{marginal cost per level 2 unit.}
  \item{cK}{marginal cost per level 3 unit.}
  \item{cost}{total cost or budget.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{level 3 sample size.}
  \item{power}{statistical power (1 - type II error).}
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing.}
  \item{nJK0}{vector with a length of three to specifiy starting values for level 1, level 2, and level 3 sample sizes.}
  \item{ncase}{number of cases to show in the output.}
  \item{gm}{grid multiplier to increase the range of sample size search for each level.}
  \item{constrain}{one of the followings can be constrained at a specified cost or value: \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{optimizer}{algorithm to find optimal sample sizes given total cost, power, or MDES. Available algorithms: \code{"auglag_cobyla"}, \code{"auglag_lbfgs"}, \code{"auglag_mma"}, or \code{"auglag_slsqp"}.}
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{P}{average proportion of level 2 units randomly assigned to treatment within level 3 units.}
  \item{g2}{number of covariates at level 2.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{R22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
}


\details{

Constrained optimal sample allocation (COSA; Hedges & Borenstein, 2014; Raudenbush, 1997; Raudenbush & Liu, 2000)
  problems are solved using \code{nloptr} (Ypma, 2014) package, an implementation of NLopt (Johnson, n.d.) in \R (R Core Team, 2016).
  More specifically, Augmented Lagrangian method is used for global optimization (AUGLAG, Birgin & Martines, 2008; Conn, Gould, & Toint, 1991)
  in conjuction with one of the following local optimization algorithms:
  Constrained Optimization by Linear Approximations (COBYLA, Powell, 1994),
  Low Storage BFGS (LBFGS, Liu & Nocedal, 1989),
  Method of Moving Asymptotes (MMA, Svanberg, 2002),
  or Sequental Least-Squares Quadratic Programming (SLSQP, Kraft, 1988).
  See \href{http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms}{NLopt website} for a brief summary.

\code{nloptr} returns values that are not integer. Rounding may produce cost, power or MDES values different from what was specified.
  A better solution is approximated using brute force. If the constrained value (cost, power or MDES)
  in the output deviates from what was specified, increasing grid multiplier (\code{gm}) often solves the problem.
  More cases can be printed by increasing \code{ncase}.

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

     optimal.bcra3f2(cn=1, cJ=10, cK=100, cost=5600,
                    constrain="cost",
                    rho2=.10)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}

\references{
Birgin, E. G., Martinez, J. M. (2008) Improving ultimate convergence of an augmented Lagrangian method. \emph{Optimization Methods and Software 23(2)}, 177-195.


Conn, A. R., Gould, N. I. M., & Toint, P.L. (1991). A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds. \emph{SIAM J. Numer. Anal. 28(2)}, 545-572.

Dong & Maynard (2013). PowerUp!: A Tool for Calculating Minum Detectable Effect Sizes and Minimum Required Sample Sizes for Experimental and Quasi-Experimental Design Studies,\emph{Journal of Research on Educational Effectiveness, 6(1)}, 24-6.

Hedges, L. V., & Borenstein, M. (2014). Conditional Optimal Design in Three- and Four-Level Experiments. \emph{Journal of Educational and Behavioral Statistics, 39(4)}, 257-281

Hedges, L. & Rhoads, C.(2009). Statistical Power Analysis in Education
Research (NCSER 2010-3006). Washington, DC: National Center for Special Education
Research, Institute of Education Sciences, U.S. Department of Education. This report is
available on the IES website at \url{http://ies.ed.gov/ncser/}.

Johnson, S. G. (n.d.). The NLopt nonlinear-optimization package. Package available at \url{http://ab-initio.mit.edu/nlopt}.

Kraft, D. (1988). A software package for sequential quadratic programming. Obersfaffeuhofen, Germany: DFVLR.

Liu, D. C., & Nocedal, J. (1989). On the limited memory BFGS method for large scale optimization. \emph{Mathematical programming, 45(1-3)}, 503-528.

Powell, M. J. (1994). A direct search optimization method that models the objective and constraint functions by linear interpolation. In Advances in optimization and numerical analysis (pp. 51-67). Springer Netherlands.

R Core Team (2016). R: A language and environment for statistical computin . R Foundation for Statistical Computing, Vienna, Austria. \url{https://www.R-project.org}.

Raudenbush, S. W. (1997). Statistical analysis and optimal
design for cluster randomized trials. \emph{Psychological Methods,
2}, 173-185.

Raudenbush, S. W., & Liu, X. (2000). Statistical power and
optimal design for multisite trials. \emph{Psychological Methods, 5}, 199-213.

Svanberg, K. (2002). A class of globally convergent optimization methods based on conservative convex separable approximations. \emph{SIAM journal on optimization, 12(2)}, 555-573.

Ypma, J. (2014). \code{nloptr}: R interface to NLopt. R package version 1.0.4. Package available at \url{https://cran.r-project.org/package=nloptr}

}
\seealso{
  \code{\link{mdes.bcra3f2}, \link{power.bcra3f2}, \link{mrss.bcra3f2}}
}
\keyword{optimal}
