
\name{power.bira2c1}
\alias{power.bira2c1}
\title{Model 2.1: Statistical Power Calculator for 2-Level Constant Effects Blocked Individual Random Assignment Designs, Treatment at Level 1}

\description{
  \code{power.bira2c1} calculates statistical power for designs with 2-levels
    where level 1 units are randomly assigned to treatment and control groups within level 2 units (school intercepts only).
}

\usage{
  power.bira2c1(mdes=.25, alpha=.05, two.tail=TRUE,
               P=.50, g1=0, R12=0,
               n, J, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{P}{average proportion of level 1 units randomly assigned to treatment within level 2 units.}
  \item{g1}{number of covariates at level 1.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{level 2 sample size.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}


\details{

Power formula was derived within power analysis framework descibed by Hedges & Rhoads (2009). Further definition of design parameters can be found in Dong & Maynard (2013).

}

\value{
  \item{fun}{function name.}
  \item{par}{list of parameters used in power calculation.}
  \item{df}{model degrees of freedom}
  \item{M}{multiplier for MDES calcuation given model degrees of freedom, \eqn{\alpha} and \eqn{\beta} (1-power).}
  \item{power}{statistical power (1 - type II error).}
}


\examples{
\dontrun{

   power.bira2c1(n=55, J=3)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}

\references{

Dong & Maynard (2013). PowerUp!: A Tool for Calculating Minum Detectable Effect Sizes and Minimum Required Sample Sizes for Experimental and Quasi-Experimental Design Studies,\emph{Journal of Research on Educational Effectiveness, 6(1)}, 24-6.

Hedges, L. & Rhoads, C.(2009). Statistical Power Analysis in Education
Research (NCSER 2010-3006). Washington, DC: National Center for Special Education
Research, Institute of Education Sciences, U.S. Department of Education. This report is
available on the IES website at \url{http://ies.ed.gov/ncser/}.

}
\seealso{
  \code{\link{mdes.bira2c1}, \link{mrss.bira2c1}, \link{optimal.bira2c1}}
}
\keyword{power}
