mdes.bcra4f3 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                        rho2, rho3, p=.50, r21=0, r22=0, r23=0, g3=0,
                        n, J, K, L){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- L * (K - 2) - g3
  SSE <- sqrt(rho3*(1-r23)/(p*(1-p)*K*L) +
                  rho2*(1-r22)/(p*(1-p)*J*K*L) +
                  (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*K*L*n))

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.bcra4f3",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3,
                                p=p, r21=r21, r22=r22, r23=r23, g3=g3,
                                n=n, J=J, K=K, L=L),
                   df=df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("main", "mdes")
  return(invisible(mdes.out))
}

# example
# mdes.bcra4f3(alpha=.05, two.tailed=TRUE, power=.80, rho3=.15, rho2=.15, n=10, J=4, K=4, L=15)

power.bcra4f3 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                         rho2, rho3, p=.50, r21=0, r22=0, r23=0, g3=0,
                         n, J, K, L){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- L * (K - 2) - g3
  SSE <- sqrt(rho3*(1-r23)/(p*(1-p)*K*L) +
                     rho2*(1-r22)/(p*(1-p)*J*K*L) +
                     (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*K*L*n))

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)
  power.out <-  list(fun = "power.bcra4f3",
                     parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                  rho2=rho2, rho3=rho3,
                                  p=p, r21=r21, r22=r22, r23=r23, g3=g3,
                                  n=n, J=J, K=K, L=L),
                     df=df,
                     ncp = es/SSE,
                     power = power)
  class(power.out) <- c("main", "power")
  return(invisible(power.out))
}
# example
# power.bcra4f3(es=0.24, alpha=.05, two.tailed=TRUE, rho3=.15, rho2=.15, n=10, J=4, K=4, L=15)

mrss.bcra4f3 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                         n, J, K, L0=10, tol=.10,
                         rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <- L0*(K-2)-g3
    if(df<=0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    L1 <- (M/es)^2 * (rho3*(1-r23)/(p*(1-p)*K) +
                          rho2*(1-r22)/(p*(1-p)*K*J) +
                          (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*K*n))
    if(abs(L1-L0)<tol){conv <- TRUE}
    L0 <- (L1+L0)/2
    i <- i+1
  }
  L <- ifelse(df>0,round(L0),NA)

  mrss.out <-  list(fun = "mrss.bcra4f3",
                    parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                 n=n, J=J, K=K, L0=L0, tol=tol,
                                 rho2=rho2, rho3=rho3,
                                 p=p, r21=r21, r22=r22, r23=r23, g3=g3),
                    df=df,
                    ncp = M,
                    L = L)
  class(mrss.out) <- c("main", "mrss")
  cat("L =", L, "\n")
  return(invisible(mrss.out))
}
# example
# mrss.bcra4f3(rho3=.15, rho2=.15, g3=1, p=.5, r23=.5, r22=.5, r21=.5, n=10, J=4, K=4)
