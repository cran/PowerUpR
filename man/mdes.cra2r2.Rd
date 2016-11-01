
\name{mdes.cra2r2}
\alias{mdes.cra2r2}
\title{Model 3.1: MDES Calculator for 2-Level Cluster Random Assignment Design, Treatment at Level 2}

\description{
  \code{mdes.cra2r2} calculates minimum detectable effect size (MDES) for designs with 2-levels
    where level 2 units are randomly assigned to treatment and control groups.
}

\usage{
  mdes.cra2r2(power=.80, alpha=.05, two.tail=TRUE,
              rho2,
              P=.50,  g2=0, R12=0, R22=0,
              n, J, ...)
}

\arguments{
  \item{power}{statistical power (1 - type II error).}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{P}{proportion of level 2 units randomly assigned to treatment.}
  \item{g2}{number of covariates at level 2.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{R22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{level 2 sample size.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}

\value{
  \item{fun}{function name.}
  \item{par}{list of parameters used in MDES calculation.}
  \item{df}{model degrees of freedom}
  \item{M}{multiplier for MDES calcuation given model degrees of freedom, \eqn{\alpha} and \eqn{\beta} (1-power).}
  \item{mdes}{minimum detectable effect size and 95\% lower and upper confidence limits.}
}

\details{

MDES formula and further definition of design parameters can be found in Dong & Maynard (2013).

}

\examples{
\dontrun{

 mdes.cra2r2(rho2=.20,
             n=4, J=20)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}

\references{

Dong & Maynard (2013). PowerUp!: A Tool for Calculating Minum Detectable Effect Sizes and Minimum Required Sample Sizes for Experimental and Quasi-Experimental Design Studies,\emph{Journal of Research on Educational Effectiveness, 6(1)}, 24-6.

}
\seealso{
  \code{\link{power.cra2r2}, \link{mrss.cra2r2}, \link{optimal.cra2r2}}
}
\keyword{mdes}
