# type I and type II error plot
.t1t2 <- function(ncp, df, alpha, two.tailed){

  # t-critical, power, beta
  talpha <- ifelse(two.tailed == FALSE,
                   qt(alpha, df, lower.tail = FALSE),
                   qt(alpha/2, df, lower.tail = FALSE)
  )
  power <- ifelse(two.tailed == FALSE,
                  1 - pt(talpha, df, ncp),
                  1 - pt(talpha, df, ncp) +
                    pt(-talpha, df, ncp)
  )
  beta = 1 - power

  # define functions
  funt0 <- function(x){
    dt(x, df = df, ncp = 0)
  } # central t
  funt1 <- function(x){
    dt(x, df = df, ncp = ncp)
  } #non-central t

  # plot central t distribution
  plot(funt0, xlim = c(-3,8), ylim = c(0, 0.5),
       yaxs = "i", xaxs = "i", bty = "l",
       sub = paste("Type I Error Rate = ", round(alpha ,digits = 3), ",",
                   "Type II Error Rate = ", round(beta, digits = 3), "\n",
                   "Non-centrality Parameter (NCP) =  ", round(ncp, digits = 3)),
       xlab = "", ylab = "")

  par(new = TRUE)

  # plot non-central t distribution
  plot(funt1, xlim = c(-3,8), ylim = c(0, 0.5),
       yaxs = "i", xaxs = "i", yaxt = "n", xaxt = "n",
       bty = "l", xlab = "", ylab = "")
  legend("topright",
         c(expression(alpha),
           expression(beta),
           ifelse(two.tailed == TRUE,
                  expression(t[alpha/2]),
                  expression(t[alpha])
           ),
           expression(NCP)),
         pch = c(19, 19, NA, NA), lty = c(NA, NA, 2, 2), cex = 1,
         col = c(adjustcolor(2, alpha.f = 0.3), adjustcolor(4, alpha.f = 0.3), 2, 4))

  # axes labels and subtitle
  title(ylab = "P(t)", line = 2)
  title(xlab = "t", line = 2)

  # draw vertical lines
  abline(v = 0, lty = 2, col = 4) # mean of central t in dashed blue line
  abline(v = ncp, lty = 2, col = 4) # mean of non-central t in dashed blue line
  abline(v = talpha, lty = 2, col = 2) # t-critical in dashed red line

  # shaded area in red for alpha
  xalpha0 <- seq(from = talpha, to = 10, by = .001)
  yalpha0 <- rep(NA, length(xalpha0))
  for(i in 1:length(xalpha0)){
    yalpha0[i] <- funt0(xalpha0[i])
  }

  xalpha <- c(xalpha0, rev(xalpha0))
  yalpha <- c(yalpha0, rep(0, length(yalpha0)))
  polygon(x = xalpha, y = yalpha, col = adjustcolor(2, alpha.f = 0.3))

  # shaded area in light blue for beta
  xbeta0 <- seq(from = -3,to = talpha, by = .001)
  ybeta0 <- rep(NA, length(xbeta0))
  for(i in 1:length(xbeta0)){
    ybeta0[i] <- funt1(xbeta0[i])
  }
  xbeta <- c(xbeta0, rev(xbeta0))
  ybeta <- c(ybeta0, rep(0, length(ybeta0)))
  polygon(x = xbeta, y = ybeta, col = adjustcolor(4, alpha.f = 0.3))
}


# Type I  and Type II error plot wrapper
t1t2.error <- function(object){

  if(inherits(object, "mrss")){
    object <- mrss.to.power(object)
  }

  # plot type I and type II error rates
  return(.t1t2(ncp = object$ncp,
               df = object$df,
               alpha = object$parms$alpha,
               two.tailed = object$parms$two.tailed))
}

# examples
# design1 <- mdes.bcra4f3(rho3=.15, rho2=.15, n=10, J=4, L=27, K=4)
# design2 <- power.bcra4f3(rho3=.15, rho2=.15, n=10, J=4, L=27, K=4)
# t1t2.error(design1)
# t1t2.error(design2)

