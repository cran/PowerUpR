# minimum detectable effect size
.mdes.fun <- function(power, alpha, sse, df, two.tailed){
  if(length(sse) > 1 || !is.numeric(sse) || sse < 0 ||
     length(df) > 1 || !is.numeric(df) || df < 1) {
    stop("Design not feasible", call. = FALSE)
  }
  t1 <- ifelse(two.tailed == TRUE, abs(qt(alpha / 2, df)), abs(qt(alpha, df)))
  t2 <- abs(qt(power, df))
  m <- ifelse(power >= 0.5, t1 + t2, t1 - t2)
  mdes <- m * sse
  lcl <- mdes * (1 - t1 / m)
  ucl <- mdes * (1 + t1 / m)
  mlu <- cbind(mdes, lcl, ucl)
  colnames(mlu) <- c("mdes", paste(100 * (1 - round(alpha, 2)), "% lcl", sep = ""),
                     paste(100 * (1 - round(alpha, 2)), "% ucl", sep = ""))
  return(invisible(mlu))
}

# statistical power
.power.fun <- function(es, alpha, sse, df, two.tailed){
  if(length(sse) > 1 || !is.numeric(sse) || sse < 0 ||
     length(df) > 1 || !is.numeric(df) || df < 1) {
    stop("Design not feasible", call. = FALSE)
  }
  lambda <- es/sse
  power <- ifelse(two.tailed == FALSE,
                  1 - pt(qt(alpha, df, lower.tail = FALSE), df, lambda),
                  1 - pt(qt(alpha / 2, df, lower.tail = FALSE), df, lambda) +
                    pt(-qt(alpha / 2, df, lower.tail = FALSE), df, lambda))
  return(invisible(power))
}


# summarize mdes output
.summ.mdes <- function(power, alpha, sse, df, two.tailed, mdes) {
  cat("\nMinimum detectable effect size: \n--------------------------------------- \nMDES is ",
      round(mdes[1], 3), " ", 100 * (1 - round(alpha, 2)), "% CI [", round(mdes[2], 3),
      ",", round(mdes[3], 3), "] with ", round(power,3)*100,
      "% power \n---------------------------------------\nDegrees of freedom: ", df,
      "\nStandardized standard error: ", round(sse, 3), "\nType I error rate: ", alpha,
      "\nType II error rate: ", round(1 - power, 3), "\nTwo-tailed test: ", two.tailed, "\n",
      sep = "")
}

# summarize power output
.summ.power <- function(es, alpha, sse, df, two.tailed, power) {
  cat("\nStatistical power: \n--------------------------------------- \n ",
      round(power, 3) * 100, "% power \n--------------------------------------- \nDegrees of freedom: ", df,
      "\nStandardized standard error: ", round(sse, 3), "\nType I error rate: ", alpha,
      "\nType II error rate: ", round(1 - power, 3), "\nTwo-tailed test: ", two.tailed, "\n",
      sep = "")
}

# path error variance functions
.se.a221 <- function(esa, r2m2, p, J) {
  var.a221 <- (1-(r2m2 + p * (1 - p) * esa^2)) / (p * (1 - p) * J)
  if(var.a221 <= 0) {
    stop("Design not feasible", call. = FALSE)
  }
  return(sqrt(var.a221))
}

.se.b221 <- function(esa, esb, escp, rho2, r22, r21, r2m2, p, n, J) {
  var.b221 <- (rho2 * (1 - (r22 + p * (1 - p) * (esa * esb + escp)^2 / rho2 + (esb^2 / rho2) * (1 - r2m2 - p * (1 - p) * esa^2))) +
            (1 - rho2) * (1 - r21) / n) / (J * (1 - (r2m2 + p * (1 - p) * esa^2)))
  if(var.b221 <= 0) {
    stop("Design not feasible", call. = FALSE)
  }
  sqrt(var.b221)
}

.se.a211 <- function(esa, rhom2, r2m1, r2m2, n, J, p) {
  t2mbar <- rhom2 * (1 - r2m2 - (p * (1 - p)  * esa^2) / rhom2)
  sig2mbar <- (1 - rhom2) * (1 - r2m1)
  var.a211 <- (t2mbar + sig2mbar / n) / (J * p * (1 - p))
  if(var.a211 <= 0) {
    stop("Design not feasible", call. = FALSE)
  }
  return(sqrt(var.a211))
}

.se.b1211 <- function(esb1, rho2, rhom2, r21, r2m1, n, J) {
  sig2mbar <- (1 - rhom2) * (1 - r2m1)
  sig2ybar <- (1 - rho2) * (1 - r21 - (((1 - rhom2) / (1 - rho2)) * esb1^2 * (1 - r2m1)))
  var.b1211 <- sig2ybar / ((J * n - J) * sig2mbar)
  if(var.b1211 <= 0) {
    stop("Design not feasible", call. = FALSE)
  }
  return(sqrt(var.b1211))
}

.se.B211 <- function(esa, esB, esb1, escp, rho2, rhom2, r22, r21, r2m2, r2m1, n, J, p) {
  t2mbar <- rhom2 * (1 - r2m2 - (p * (1 - p) * esa^2) / rhom2)
  sig2mbar <- (1 - rhom2) * (1 - r2m1)
  t2ybar <- rho2 * (1 - r22) - p * (1 - p) * (esa * esB + escp)^2 -
    ((1 / (p * (1 - p))) * esB^2 * rhom2 * (1 - r2m2) +
       (1 / (p * (1 - p))) * esB^2 * (1 - rhom2) * (1 - r2m1) / n - esa^2 * esB^2) / (1 / (p * (1 - p)))
  sig2ybar <- (1 - rho2) * (1 - r21 - (((1 - rhom2) / (1 - rho2)) * esb1^2 * (1 - r2m1)))
  var.B211 <- (t2ybar + sig2ybar / n) / (J * (t2mbar + sig2mbar / n))
  if(var.B211 <= 0) {
    stop("Design not feasible", call. = FALSE)
  }
  return(sqrt(var.B211))
}

.se.sobel <- function(x, y, sex, sey) {
  var.sobel <- (x^2 * sey^2  + y^2 * sex^2)
  if(var.sobel <= 0) {
    stop("Design not feasible", call. = FALSE)
  }
  return(sqrt(var.sobel))
}

# power functions
.power.sobel <- function(x, y, sex, sey, alpha, two.tailed, df = 1e+8) {
  sesobel <- .se.sobel(x = x, y = y, sex = sex, sey = sey)
  capture.output({
    power <- .power.fun(es = x*y, alpha = alpha, sse = sesobel, two.tailed = two.tailed, df = df)
  })
  return(power)
}

.power.jt <- function(x, y, sex, sey, alpha, two.tailed, dfx, dfy) {
  capture.output({
    powerx <- .power.fun(es = x, alpha = alpha, sse = sex, two.tailed = two.tailed, df = dfx)
    powery <- .power.fun(es = y, alpha = alpha, sse = sey, two.tailed = two.tailed, df = dfy)
  })
  return( powerx*powery)
}

.power.mc <- function(nsims, ndraws, x, y, sex, sey, alpha, two.tailed) {
  rejmc <- NULL
  for (i in 1:nsims){
    xstar <- rnorm(1, x, sex)
    ystar <- rnorm(1, y, sey)
    rejmc <- c(rejmc, quantile(rnorm(ndraws, xstar, sex)*rnorm(ndraws, ystar, sey), probs = ifelse(two.tailed, alpha/2, alpha), na.rm = TRUE) > 0)
  }
  return(mean(rejmc))
}

