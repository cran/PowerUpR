# type I and type II error plot
.t1t2 <- function(ncp, df, alpha, two.tail){

  # t-critical, power, beta
  talpha <- ifelse(two.tail == FALSE,
                   qt(alpha, df, lower.tail = FALSE),
                   qt(alpha/2, df, lower.tail = FALSE)
  )
  power <- ifelse(two.tail == FALSE,
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
  plot(funt0, xlim = c(-3,10), ylim = c(0, 0.5),
       yaxs = "i", xaxs = "i", bty = "l",
       main = expression(paste("Type I ", (alpha), " and "," Type II ", (beta), " Error Rates")),
       sub = paste("Type I = ", round(alpha ,digits = 3), ",",
                   "Type II = ", round(beta, digits = 3), ",",
                   "ncp =  ", round(ncp, digits = 3)),
       xlab = "", ylab = "")

  par(new = TRUE)

  # plot non-central t distribution
  plot(funt1, xlim = c(-3,10), ylim = c(0, 0.5),
       yaxs = "i", xaxs = "i", yaxt = "n", xaxt = "n",
       bty = "l", xlab = "", ylab = "")
  legend("topright",
         c(expression(alpha),
           expression(beta),
           ifelse(two.tail == TRUE,
                  expression(t[alpha/2]),
                  expression(t[alpha])
           ),
           expression(ncp)),
         x.intersp = 0.9, y.intersp = 0.9,
         pch = c(19, 19, NA, NA), lty = c(NA, NA, 2, 2), cex = 1,
         col = c(adjustcolor(2, alpha.f = 0.3), adjustcolor(4, alpha.f = 0.3), 2, 4),
         bty = 1)

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
.t1t2.error <- function(x){

  if(class(x)[2] == "cosa"){
    x <- cosa.to.power(x)
  }

  # plot type I and type II error rates
  return(.t1t2(ncp = x$ncp,
               df = x$parms$df,
               alpha = x$parms$alpha,
               two.tail = x$parms$two.tail))
}
