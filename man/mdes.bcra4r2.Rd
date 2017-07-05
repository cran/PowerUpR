\name{mdes.bcra4r2}
\alias{mdes.bcra4r2}
\title{MDES Calculator for Four-Level Random Effects Block Random Assignment Designs, Treatment at Level-2}

\description{
  \code{mdes.bcra4r2} calculates minimum detectable effect size (MDES) for designs with four-levels
    where level-2 units are randomly assigned to treatment and control groups within level-3 units (random blocks).
}

\usage{
  mdes.bcra4r2(power=.80, alpha=.05, two.tail=TRUE,
       rho2, rho3, rho4, omega3, omega4,
       P=.50, R12=0, R22=0, RT32=0, RT42=0, g4=0,
       n, J, K, L, ...)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{rho4}{proportion of variance in the outcome explained by level-4 units.}
  \item{omega3}{treatment effect heterogeneity as ratio of treatment effect variance among level-3 units to the residual variance at level-3.}
  \item{omega4}{treatment effect heterogeneity as ratio of treatment effect variance among level-4 units to the residual variance at level-4.}
  \item{P}{average proportion of level-2 units randomly assigned to treatment within level-3 units.}
  \item{g4}{number of covariates at level-4.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{RT32}{proportion of treatment effect variance among level-3 units explained by level-3 covariates.}
  \item{RT42}{proportion of treatment effect variance among level-4 units explained by level-4 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{harmonic mean of level-2 units across level-3 units (or simple average).}
  \item{K}{harmonic mean of level-3 units across level-4 units (or simple average).}
  \item{L}{number of level-4 units.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in MDES calculation.}
  \item{ncp}{noncentrality parameter.}
  \item{mdes}{minimum detectable effect size and \eqn{(1-\alpha)}\% lower and upper confidence limits.}
}

\examples{
\dontrun{

  mdes.bcra4r2(rho4=.05, rho3=.15, rho2=.15,
              omega4=.50, omega3=.50,
              n=10, J=4, L=27, K=4)
  }
}

\seealso{
  \code{\link{power.bcra4r2}, \link{cosa.bcra4r2}}
}

\keyword{mdes}
\keyword{main effect}
\keyword{cluster random assignment}
\keyword{block design}