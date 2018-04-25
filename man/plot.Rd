\title{Plots}
\name{plots}
\alias{plot.power}
\alias{plot.mdes}

\description{
  \code{plot} plots power or MDES curves.
}

\usage{
  \method{plot}{power}(x, ypar = "power", xpar, ylim, xlim, ...)
  \method{plot}{mdes}(x, ypar = "power", xpar, ylim, xlim, ...)
}

\arguments{
  \item{x}{an object returned from one of the \code{PowerUpR} functions.}
  \item{ypar}{character; \code{"mdes"}, or \code{"power"}.}
  \item{xpar}{character; one of the sample sizes.}
  \item{ylim}{limits for \code{ypar}.}
  \item{xlim}{limits for \code{xpar}.}
  \item{...}{other graphical parameters.}
}


\examples{
design1 <- mdes.cra3r3(rho3=.06, rho2=.17, n=15, J=3, K=60)
plot(design1, ypar = "mdes", xpar = "K", xlim = c(30, 100))
plot(design1, ypar = "power", xpar = "K", xlim = c(30, 100))

design2 <- power.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3, K=60)
plot(design2, ypar = "mdes", xpar = "K", xlim = c(30, 100))
plot(design2, ypar = "power", xpar = "K", xlim = c(30, 100))
}
