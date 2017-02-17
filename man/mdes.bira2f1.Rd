
\name{mdes.bira2f1}
\alias{mdes.bira2f1}
\title{Model 2.2: MDES Calculator for 2-Level Fixed Effects Blocked Individual Random Assignment Designs, Treatment at Level 1}

\description{
  \code{mdes.bira2f1} calculates minimum detectable effect size (MDES) for designs with 2-levels
    where level 1 units are randomly assigned to treatment and control groups within level 2 units (fixed blocks).
}

\usage{
  mdes.bira2f1(power=.80, alpha=.05, two.tail=TRUE,
              P=.50, g1=0, R12=0,
              n, J, ...)
}

\arguments{
  \item{power}{statistical power (1 - type II error).}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{P}{average proportion of level 1 units randomly assigned to treatment within level 2 units.}
  \item{g1}{number of covariates at level 1.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{level 2 sample size.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}

\value{
  \item{fun}{function name.}
  \item{par}{list of parameters used in MDES calculation.}
  \item{df}{degrees of freedom}
  \item{M}{multiplier for MDES calcuation given degrees of freedom, \eqn{\alpha} and \eqn{\beta} (1-power).}
  \item{mdes}{minimum detectable effect size and 95\% lower and upper confidence limits.}
}

\details{

MDES formula and further definition of design parameters can be found in Dong & Maynard (2013).

}

\examples{
\dontrun{

 mdes.bira2f1(n=55, J=3)

  }
}

\author{
  Metin Bulus \email{bulus.metin@gmail.com}
  Nianbo Dong \email{dong.nianbo@gmail.com}
}

\references{

Dong, N., & Maynard, R. A. (2013). PowerUp!: A Tool for Calculating Minum Detectable Effect Sizes and Minimum Required Sample Sizes for Experimental and Quasi-Experimental Design Studies,\emph{Journal of Research on Educational Effectiveness, 6(1)}, 24-6.

}
\seealso{
  \code{\link{power.bira2f1}, \link{mrss.bira2f1}, \link{optimal.bira2f1}}
}
\keyword{mdes}
