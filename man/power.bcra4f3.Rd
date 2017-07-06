
\name{power.bcra4f3}
\alias{power.bcra4f3}
\title{Statistical Power Calculator for Four-Level Fixed Effects Blocked Cluster Random Assignment Designs, Treatment at Level-3}

\description{
  \code{power.bcra4f3} calculates statistical power for designs with four-levels
    where level-3 units are randomly assigned to treatment and control groups within level-4 units (fixed blocks).
}

\usage{
  power.bcra4f3(mdes=.25, alpha=.05, two.tail=TRUE,
       rho2, rho3,
       P=.50, R12=0, R22=0, R32=0, g3=0,
       n, J, K, L, ...)
}

\arguments{
  \item{mdes}{minimum detectable effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tail}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome explained by level-2 units.}
  \item{rho3}{proportion of variance in the outcome explained by level-3 units.}
  \item{P}{average proportion of level-3 units randomly assigned to treatment within level-4 units.}
  \item{g3}{number of covariates at level-3.}
  \item{R12}{proportion of level-1 variance in the outcome explained by level-1 covariates.}
  \item{R22}{proportion of level-2 variance in the outcome explained by level-2 covariates.}
  \item{R32}{proportion of level-3 variance in the outcome explained by level-3 covariates.}
  \item{n}{harmonic mean of level-1 units across level-2 units (or simple average).}
  \item{J}{harmonic mean of level-2 units across level-3 units (or simple average).}
  \item{K}{harmonic mean of level-3 units across level-4 units (or simple average).}
  \item{L}{number of level-4 units.}
  \item{...}{to handle extra parameters passed from other functions, do not define any additional parameters.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
}


\examples{
\dontrun{

    power.bcra4f3(rho3=.15, rho2=.15,
                 n=10, J=4, K=4, L=15)
				 
  }
}

\seealso{
  \code{\link{mdes.bcra4f3}, \link{optimal.bcra4f3}}
}

\keyword{power}
\keyword{main effect}
\keyword{cluster random assignment}
\keyword{block design}
