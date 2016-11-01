
\name{power.bcra3f2}
\alias{power.bcra3f2}
\title{Model 4.1: Statistical Power Calculator for 3-Level Fixed Effects Blocked Cluster Random Assignment Design, Treatment at Level 2}

\description{
  \code{power.bcra3f2} calculates statistical power for designs with 3-levels
    where level 2 units are randomly assigned to treatment and control groups within level 3 units (fixed blocks).
}

\usage{
  power.bcra3f2(mdes=.25, alpha=.05, two.tail=TRUE,
       rho2,
       P=.50, g2=0, R12=0, R22=0,
       n, J, K, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{P}{average proportion of level 2 units randomly assigned to treatment within level 3 units.}
  \item{g2}{number of covariates at level 2.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{R22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{level 3 sample size.}
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

    power.bcra3f2(rho2=.10,
                 n=20, J=44, K=5)

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
  \code{\link{mdes.bcra3f2}, \link{mrss.bcra3f2}, \link{optimal.bcra3f2}}
}
\keyword{power}
