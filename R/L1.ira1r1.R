
mdes.ira1r1 <- function(power=.80, alpha=.05, two.tail=TRUE,
                             P=.50, g1=0, R12=0, n, ...){
  df <- n-g1-2
  SSE <- sqrt((1-R12)/(P*(1-P)*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.ira1r1",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}
# example
# mdes.ira1r1(n=200)


power.ira1r1 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                              P=.50, g1=0, R12=0, n, ...){
  df <- n-g1-2
  SSE <- sqrt((1-R12)/(P*(1-P)*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.ira1r1",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.ira1r1(n=200)

cosa.ira1r1 <- function(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
                        n0=10, tol=.10,
                        P=.50, g1=0, R12=0){
    parms <- as.list(environment(), all=TRUE)
    .error.handler(parms)
    i <- 0
    conv <- FALSE
    while(i<=100 & conv==FALSE){
      df <- n0-g1-2
      if(df<= 0 | is.infinite(df)){break}
      T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
      T2 <- abs(qt(power,df))
      M <- ifelse(power>=.5,T1+T2,T1-T2)
      n1 <- (M/mdes)^2 * ((1-R12)/(P*(1-P)))
      if(abs(n1-n0)<tol){conv <- TRUE}
      n0 <- (n1+n0)/2
      i <- i+1
    }
    n <- round(ifelse(df>0,round(n0),NA))
    n.out <-  list(fun = "cosa.ira1r1",
                       parms = parms,
                       n = n)
    class(n.out) <- c("parms", "cosa")

  print(paste("n =",round(n)))
  return(invisible(n.out))
}

optimal.ira1r1 <- function(...){
  .Deprecated(new="cosa.ira1r1")
  cosa.ira1r1(...)
}

# example
# cosa.ira1r1()
