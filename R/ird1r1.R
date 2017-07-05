mdes.ird1r1 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                        g1=0, R12=0, n, ...){
  df <- n - g1 - 2
  D <- .D.fun(P, k1, k2, dist.Z, RTZ)
  SSE <- sqrt(D*(1-R12)/(P*(1-P)*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.ird1r1",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.ird1r1(n=400)

power.ird1r1 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                        P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                        g1=0, R12=0,  n, ...){
  df <- n - g1 - 2
  D <- .D.fun(P, k1, k2, dist.Z, RTZ)
  SSE <- sqrt(D*(1-R12)/(P*(1-P)*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.ird1r1",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.ird1r1(mdes=0.466, n=400)

cosa.ird1r1 <- function(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
                        n0=10, tol=.10, P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                        g1=0, R12=0){
   D <- .D.fun(P, k1, k2, dist.Z, RTZ)
   parms <- as.list(environment(), all=TRUE)
    .error.handler(parms)
    i <- 0
    conv <- FALSE
    while(i<=100 & conv==FALSE){
      df <- n0 - g1 - 2
      if(df<= 0 | is.infinite(df)){break}
      T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
      T2 <- abs(qt(power,df))
      M <- ifelse(power>=.5,T1+T2,T1-T2)
      n1 <- (M/mdes)^2 * (D*(1-R12)/(P*(1-P)*n0))
      if(abs(n1-n0)<tol){conv <- TRUE}
      n0 <- (n1+n0)/2
      i <- i+1
    }
    n <- round(ifelse(df>0,round(n0),NA))
    n.out <-  list(fun = "cosa.ird1r1",
                       parms = parms,
                       n = n)
    class(n.out) <- c("parms", "cosa")

  print(paste("n =",round(n)))
  return(invisible(n.out))
}

# examples
# cosa.ird1r1(mdes=.466)




