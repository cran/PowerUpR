\name{cra3r3}
\alias{cra3r3}
\alias{mdes.cra3r3}
\alias{mdes.mod331}
\alias{mdes.mod332}
\alias{mdes.mod333}
\alias{mdesd.mod331}
\alias{mdesd.mod332}
\alias{mdesd.mod333}
\alias{power.cra3r3}
\alias{power.mod331}
\alias{power.mod332}
\alias{power.mod333}
\alias{mrss.cra3r3}
\alias{mrss.mod331}
\alias{mrss.mod332}
\alias{mrss.mod333}
\title{Three-level Cluster-randomized Trials to Detect Main and Moderation Effects}

\description{
   Use \code{mdes.<design>()} to calculate the minimum detectable effect size for the main effect, \code{mdesd.<design>()} to calculate the minimum detectable effect size difference for the moderation effect, \code{power.<design>()} to calculate the statistical power, and \code{mrss.<design>()} to calculate the minimum required sample size.
   Use \code{<output>.cra3r3()} for the main effect, \code{<output>.mod331()} for the moderator effect at level 1, \code{<output>.mod332()} for the moderator effect at level 2, \code{<output>.mod333()} for the moderator effect at level 3.
}

\usage{
mdes.cra3r3(power=.80, alpha=.05, two.tailed=TRUE,
            rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0,
            n, J, K)

mdesd.mod331(power=.80, alpha=.05, two.tailed=TRUE,
            rho2, rho3, omegam2=0, omegam3=0,
            g1=0, r21=0, r2m2=0, r2m3=0,
            p=.50, q=NULL, n, J, K)

mdesd.mod332(power=.80, alpha=.05, two.tailed=TRUE,
            rho2, rho3, omegam3, g2=0, r21=0, r22=0, r2m3=0,
            p=.50, q=NULL, n, J, K)

mdesd.mod333(power=.80, alpha=.05, two.tailed=TRUE,
            rho2, rho3,  g3=0, r21=0, r22=0, r23=0,
            p=.50, q=NULL, n, J, K)

power.cra3r3(es=.25, alpha=.05, two.tailed=TRUE,
             rho2, rho3, g3=0, r21=0, r22=0, r23=0,
             p=.50, q=NULL, n, J, K)

power.mod331(es=.25, alpha=.05, two.tailed=TRUE,
             rho2, rho3, omegam2, omegam3,
             g1=0, r21=0, r2m2=0, r2m3=0,
             p=.50, q=NULL, n, J, K)

power.mod332(es=.25, alpha=.05, two.tailed=TRUE,
             rho2, rho3, omegam3, g2=0, r21=0, r22=0, r2m3=0,
             p=.50, q=NULL, n, J, K)

power.mod333(es=.25, alpha=.05, two.tailed=TRUE,
             rho2, rho3, g3=0, r21=0, r22=0, r23=0,
             p=.50, q=NULL, n, J, K)

mrss.cra3r3(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            n, J, K0=10, tol=.10,
            rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0)

mrss.mod331(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            rho2, rho3, omegam2, omegam3,
            g1=0, r21=0, r2m2=0, r2m3=0,
            p=.50, q=NULL, n, J, K0=10, tol=.10)

mrss.mod332(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            rho2, rho3, omegam3, g2=0, r21=0, r22=0, r2m3=0,
            p=.50, q=NULL, n, J, K0=10, tol=.10)

mrss.mod333(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
            rho2, rho3, g3=0, r21=0, r22=0, r23=0,
            p=.50, q=NULL, n, J, K0=10, tol=.10)
}

\arguments{
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{es}{effect size.}
  \item{alpha}{probability of type I error.}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing, \code{FALSE} for one-tailed hypothesis testing. }
  \item{rho2}{proportion of variance in the outcome between level 2 units (unconditional ICC2).}
  \item{rho3}{proportion of variance in the outcome between level 3 units (unconditional ICC3).}
  \item{omegam2}{ratio of the unconditional variance in the moderator effect that is between level 2 units to the residual variance between level 2 units in the null model.}
  \item{omegam3}{ratio of the unconditional variance in the moderator effect that is between level 3 units to the residual variance between level 3 units in the null model.}
  \item{p}{proportion of level 3 units randomly assigned to treatment.}
  \item{q}{proportion of level 1, level 2, or level 3 units in the moderator subgroup.}
  \item{g1}{number of covariates at level 1.}
  \item{g2}{number of covariates at level 2.}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r23}{proportion of level 3 variance in the outcome explained by level 3 covariates.}
  \item{r2m2}{proportion of variance in the moderator effect that is explained by level 2 covariates.}
  \item{r2m3}{proportion of variance in the moderator effect that is explained by level 3 covariates.}
  \item{n}{harmonic mean of level 1 units across level 2 units (or simple average).}
  \item{J}{harmonic mean of level 2 units across level 3 units (or simple average).}
  \item{K}{level 3 sample size.}
  \item{K0}{starting value for \code{K}.}
  \item{tol}{tolerance to end iterative process for finding \code{K}.}
}

\value{
  \item{fun}{function name.}
  \item{parms}{list of parameters used in power calculation.}
  \item{df}{degrees of freedom.}
  \item{ncp}{noncentrality parameter.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mdes}{minimum detectable effect size.}
  \item{K}{number of level 3 units.}
}


\examples{
# cross-checks for the main effect
mdes.cra3r3(rho3=.06, rho2=.17, n=15, J=3, K=60)
power.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3, K=60)
mrss.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3)

# cross-checks for the randomly varying cont. L1 moderator effect
mdes.mod331(power=.80, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=.10, omegam3=.10,
            q=NULL, n=15, J=3, K=60)
power.mod331(es=0.1248, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=.10, omegam3=.10,
            q=NULL, n=15, J=3, K=60)
mrss.mod331(es=0.1248, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=.10, omegam3=.10,
            q=NULL, n=15, J=3)

# cross-checks for the non-randomly varying cont. L1 moderator effect
mdesd.mod331(power=.80, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=0, omegam3=0,
            q=NULL, n=15, J=3, K=60)
power.mod331(es=.0946, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=0, omegam3=0,
            q=NULL, n=15, J=3, K=60)
mrss.mod331(es=.0946, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=0, omegam3=0,
            q=NULL, n=15, J=3)

# cross-checks for the randomly varying bin. L1 moderator effect
mdesd.mod331(power=.80, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=.10, omegam3=.10,
            q=.50, n=15, J=3, K=60)
power.mod331(es=.2082, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=.10, omegam3=.10,
            q=.50, n=15, J=3, K=60)
mrss.mod331(es=.2082, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=.10, omegam3=.10,
            q=.50, n=15, J=3)

# cross-checks for the non-randomly varying bin. L1 moderator effect
mdesd.mod331(power=.80, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=0, omegam3=0,
            q=.50, n=15, J=3, K=60)
power.mod331(es=.1893, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=0, omegam3=0,
            q=.50, n=15, J=3, K=60)
mrss.mod331(es=.1893, alpha=.05, two.tailed=TRUE,
            rho2=.17, rho3=.06, omegam2=0, omegam3=0,
            q=.50, n=15, J=3)
}

\seealso{
   For a more flexible sample size determination see \code{\link[cosa]{cosa.crd3r3}.}
}
