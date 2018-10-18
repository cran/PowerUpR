\title{Plots}
\name{plots}
\alias{plot.power}
\alias{plot.mdes}
\alias{plot.mrss}

\description{Plots statistical power or minimum detectable effect size curves with (1-\eqn{\alpha})x100 \% confidence interval.
}

\usage{
  \method{plot}{power}(x, ypar = "mdes",  xpar = NULL,
          xlim = NULL, ylim = NULL,
          xlab = NULL, ylab = NULL,
          main = NULL, sub = NULL,
          locate = FALSE, ...)

  \method{plot}{mdes}(x, ypar = "mdes",  xpar = NULL,
          xlim = NULL, ylim = NULL,
          xlab = NULL, ylab = NULL,
          main = NULL, sub = NULL,
          locate = FALSE, ...)

  \method{plot}{mrss}(x, ypar = "mdes",  xpar = NULL,
          xlim = NULL, ylim = NULL,
          xlab = NULL, ylab = NULL,
          main = NULL, sub = NULL,
          locate = FALSE, ...)
}

\arguments{
  \item{x}{an object returned from one of the \code{PowerUpR} functions.}
  \item{ypar}{character; \code{"mdes"} or  \code{"power"} on \code{y} axis .}
  \item{xpar}{character; one of the sample sizes on \code{x} axis.}
  \item{xlim}{limits for \code{xpar}.}
  \item{ylim}{limits for \code{ypar}.}
  \item{xlab}{\code{x} axis label.}
  \item{ylab}{\code{y} axis label.}
  \item{main}{title for the plot.}
  \item{sub}{subtitle for the plot.}
  \item{locate}{logical; \code{TRUE} locates parameter values for design \code{x} on the plot.}
  \item{...}{other graphical parameters to pass to \code{plot.new()}.}
}


\examples{
design1 <- mdes.cra3r3(rho3=.06, rho2=.17, n=15, J=3, K=60)
plot(design1, ypar = "mdes", xpar = "K", xlim = c(30, 100))
plot(design1, ypar = "power", xpar = "K", xlim = c(30, 100))

design2 <- power.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3, K=60)
plot(design2, ypar = "mdes", xpar = "K", xlim = c(30, 100))
plot(design2, ypar = "power", xpar = "K", xlim = c(30, 100))
}
