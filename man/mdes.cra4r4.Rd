
\name{mdes.cra4r4}
\alias{mdes.cra4r4}
\title{Model 3.3: MDES Calculator for 4-Level Cluster Random Assignment Designs, Treatment at Level 4}

\description{
  \code{mdes.cra4r4} calculates minimum detectable effect size (MDES) for designs with 4-levels
    where level 4 units are randomly assigned to treatment and control groups.
}

\usage{
  mdes.cra4r4(power=.80, alpha=.05, two.tail=TRUE,
       rho2, rho3, rho4,
       P=.50, R12=0, R22=0, R32=0, R42=0, g4=0,
       n, J, K, L, ...)
}

\arguments{
  \item{power}{statistical power (1 - type II error).}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level 2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level 3 units.}
  \item{rho4}{proportion of variance in the outcome explained by level 4 units.}
  \item{P}{proportion of level 4 units randomly assigned to treatment.}
  \item{g4}{number of covariates at level 4.}
  \item{R12}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{R22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{R32}{proportion of level 3 variance in the outcome explained by level 3 covariates.}
  \item{R42}{proportion of level 4 variance in the outcome explained by level 4 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{harmonic mean of level 3 units across level 3 units (or simple average).}
  \item{L}{number of level 4 units.}
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

  mdes.cra4r4(rho4=.05, rho3=.05, rho2=.10,
              n=10, J=2, K=3, L=20)
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
  \code{\link{power.cra4r4}, \link{mrss.cra4r4}, \link{optimal.cra4r4}}
}
\keyword{mdes}
