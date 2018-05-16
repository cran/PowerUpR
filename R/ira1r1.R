
mdes.ira1r1 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                             p=.50, g1=0, r21=0, n){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- n-g1-2
  SSE <- sqrt((1-r21)/(p*(1-p)*n))

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.ira1r1",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                p=p, r21=r21, g1=g1,
                                n=n),
                   df=df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("main", "mdes")
  return(invisible(mdes.out))
}
# example
# mdes.ira1r1(n=200)


power.ira1r1 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                              p=.50, g1=0, r21=0, n){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- n-g1-2
  SSE <- sqrt((1-r21)/(p*(1-p)*n))

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)
  power.out <-  list(fun = "power.ira1r1",
                     parms =list(es=es, alpha=alpha, two.tailed=two.tailed,
                                 p=p, r21=r21, g1=g1,
                                 n=n),
                     df=df,
                     ncp = es/SSE,
                     power = power)
  class(power.out) <- c("main", "power")
  return(invisible(power.out))
}
# example
# power.ira1r1(n=200)

mrss.ira1r1 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                        n0=10, tol=.10,
                        p=.50, g1=0, r21=0){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <- n0-g1-2
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    n1 <- (M/es)^2 * ((1-r21)/(p*(1-p)))
    if(abs(n1-n0)<tol){conv <- TRUE}
    n0 <- (n1+n0)/2
    i <- i+1
  }
  n <- round(ifelse(df>0,round(n0),NA))

  n.out <-  list(fun = "mrss.ira1r1",
                 parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                              n0=n0, tol=tol,
                              p=p, r21=r21, g1=g1),
                 df=df,
                 ncp = M,
                 n = n)
  class(n.out) <- c("main", "mrss")
  cat("n =", n, "\n")
  return(invisible(n.out))
}

# mrss.ira1r1()
