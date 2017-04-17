mdes.bcra4f3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, rho3, P=.50, R12=0, R22=0, R32=0, g3=0,
                        n, J, K, L, ...){
  df <- L*(K-2)-g3
  T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
  T2 <- abs(qt(power,df))
  M <- ifelse(power>=.5,T1+T2,T1-T2)
  mdes <- M*sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                  rho2*(1-R22)/(P*(1-P)*J*K*L) +
                  (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
  LCI <- mdes*(1 - T1/M)
  UCI <- mdes*(1 + T1/M)
  MLU <- cbind(mdes, LCI, UCI)
  colnames(MLU) <- c("mdes", "95% LCL", "95% UCL")
  fun <- "mdes.bcra4f3"
  par <- list(power=power, alpha=alpha, two.tail=two.tail,
              rho2=rho2, rho3=rho3,  
              P=P, g3=g3, R12=R12, R22=R22, R32=R32,
              n=n, J=J, K=K, L=L)
  mdes.out <- list(fun=fun,par=par,df=df,M=M,mdes=round(MLU,digits=3))
  class(mdes.out) <- c("pars")
  print(round(MLU, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.bcra4f3(alpha=.05, two.tail=TRUE, power=.80, rho3=.15, rho2=.15, n=10, J=4, K=4, L=15)

power.bcra4f3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, rho3, P=.50, R12=0, R22=0, R32=0, g3=0,
                         n, J, K, L, ...){
  df <- L*(K-2)-g3
  lamda <- mdes/sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                     rho2*(1-R22)/(P*(1-P)*J*K*L) +
                     (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
  power <- ifelse(two.tail==FALSE,
                  1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                  1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                    pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
  fun <- "power.bcra4f3"
  par <- list(mdes=mdes, alpha=alpha, two.tail=two.tail, 
              rho2=rho2, rho3=rho3,
              P=P, g3=g3, R12=R12, R22=R22, R32=R32, 
              n=n, J=J, K=K, L=L)
  power.out <- list(fun=fun,par=par,df=df,lamda=lamda,power=round(power,digits=3))
  class(power.out) <- c("pars")
  print(paste("power = ", round(power, digits=3)))
  return(invisible(power.out))
}
# example
# power.bcra4f3(mdes=0.24, alpha=.05, two.tail=TRUE, rho3=.15, rho2=.15, n=10, J=4, K=4, L=15)



# gm = grid multiplier: default 2
# ncase = numer of fixed samples in the output
# constrain = "power", or "mdes" :  default "power"
mrss.bcra4f3 <- function(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
                        gm=2, ncase=10, constrain="power",
                        n=NULL, J=NULL, K=NULL, L=NULL, L0=10, K0=10, tol=.10,
                        rho2, rho3, P=.50, g3=0, R12=0, R22=0, R32=0){
  
  # sample size for one level allowed to be null, the rest must be specified
  check.NULL <- c(is.null(n),is.null(J),is.null(K),is.null(L))
  if(length(check.NULL[check.NULL==TRUE])>1){
    stop("Specify three of the n, J, K, L.")
  }
  
  # equality constraints on power
  eq.power <- function(nJKL){
    n <- nJKL[1]
    J <- nJKL[2]
    K <- nJKL[3]
    L <- nJKL[4]
    df <- (L*(K-2)-g3)[[1]]
    lamda <- mdes/sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                       rho2*(1-R22)/(P*(1-P)*J*K*L) +
                       (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
    power1 <- ifelse(two.tail==FALSE,
                    1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                    1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                      pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
    power0 <- power1-power
    return(power0)
  }
  
  # equality constraints on mdes
  eq.mdes <- function(nJKL){
    n <- nJKL[1]
    J <- nJKL[2]
    K <- nJKL[3]
    L <- nJKL[4]
    df <- (L*(K-2)-g3)[[1]]
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    mdes1 <- M*sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                   rho2*(1-R22)/(P*(1-P)*J*K*L) +
                  (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
    mdes0 <-mdes1-mdes
    return(mdes0)
  }
  
  # find minimum required sample size per higher unit
  if(is.null(L)){
    i <- 0
    conv <- FALSE
    while(i<=100 & conv==FALSE){
      df <- L0*(K-2)-g3
      if(df<=0 | is.infinite(df)){break}
      T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
      T2 <- abs(qt(power,df))
      M <- ifelse(power>=.5,T1+T2,T1-T2)
      L1 <- (M/mdes)^2 * (rho3*(1-R32)/(P*(1-P)*K) +
                          rho2*(1-R22)/(P*(1-P)*K*J) +
                         (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n))
      if(abs(L1-L0)<tol){conv <- TRUE}
      L0 <- (L1+L0)/2
      i <- i+1
    }
    L <- ifelse(df>0,round(L0),NA)
  }else if(is.null(K)){
    i <- 0
    conv <- FALSE
    while(i<=100 & conv==FALSE){
      df <- L*(K0-2)-g3
      if(df<= 0 | is.infinite(df)){break}
      T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
      T2 <- abs(qt(power,df))
      M <- ifelse(power>=.5,T1+T2,T1-T2)
      K1 <- (rho3*(1-R32) + rho2*(1-R22)/J +
            (1-rho3-rho2)*(1-R12)/(n*J))/
            (P*(1-P)*L*(mdes/M)^2)
      if(abs(K1-K0)<tol){conv <- TRUE}
      K0 <- (K1+K0)/2
      i <- i+1
    }
    K <- ifelse(df>0,round(K0),NA)
  }else if(is.null(J)){
    df <- L*(K-2)-g3
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    J <- (rho2*(1-R22) + (1-rho3-rho2)*(1-R12)/n)/
         (P*(1-P)*K*L*(mdes/M)^2 - rho3*(1-R32))
    J <- round(J)
  }else if(is.null(n)){
    df <- L*(K-2)-g3
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    n <- (1-rho3-rho2)*(1-R12)/
         (P*(1-P)*J*K*L*(mdes/M)^2 -
          J*rho3*(1-R32) - rho2*(1-R22))
    n <- round(n)
  }
  
  nJKL1 <- c(n,J,K,L)
  if(any(nJKL1<=0)|any(is.na(nJKL1))){stop("Solution not feasible due to nonpositive sample size estimation.")}
  
  # round solution
  est.mdes <- eq.mdes(nJKL1) + mdes
  est.power <- eq.power(nJKL1) + power
  round.mrss <- cbind(nJKL1[1], nJKL1[2], nJKL1[3], nJKL1[4], est.mdes, est.power)
  colnames(round.mrss) <- c("n", "J", "K", "L", "mdes", "power")
 
  # grid to approximate an integer solution
  if(gm<=0){stop("'gm' should be >= 1")}
  if(ncase<=1){stop("'ncase' should be >= 2")}
  gridnJKL <- expand.grid(c(nJKL1[1]-gm*4):(nJKL1[1]+gm*4),
                          c(nJKL1[2]-gm*3):(nJKL1[2]+gm*3),
                          c(nJKL1[3]-gm*2):(nJKL1[3]+gm*2),
                          c(nJKL1[4]-gm*1):(nJKL1[4]+gm*1))
  
  # ensure nonzero postive cells while considering cells associated with df
  gridnJKL <- as.matrix(gridnJKL[gridnJKL[,1]>0 &
                                   gridnJKL[,2]>0 &
                                   gridnJKL[,3]>2 &
                                   gridnJKL[,4]>g3, ])
  
  # difference from specified power (default, .80), mdes (default, .25), or cost
  integer.mrss <- matrix(NA,nrow(gridnJKL),6)
  integer.mrss[,1:4] <- gridnJKL
  for(i in 1:nrow(gridnJKL)){
    integer.mrss[i,5] <- eq.mdes(gridnJKL[i,])
    integer.mrss[i,6] <- eq.power(gridnJKL[i,])
  }
  
  # output result
  if(constrain=="power"){
    idx <- order(abs(integer.mrss[,6]),
                 decreasing=FALSE)[1:ncase]
    integer.mrss <- integer.mrss[idx,]
    integer.mrss[,5] <- integer.mrss[,5] + mdes
    integer.mrss[,6] <- integer.mrss[,6] + power
    colnames(integer.mrss) <- c("n", "J", "K", "L", "mdes", "power")
  }else if(constrain=="mdes"){
    idx <- order(abs(integer.mrss[,5]),
                 decreasing=FALSE)[1:ncase]
    integer.mrss <- integer.mrss[idx,]
    integer.mrss[,5] <- integer.mrss[,5] + mdes
    integer.mrss[,6] <- integer.mrss[,6] + power
    colnames(integer.mrss) <- c("n", "J", "K", "L", "mdes", "power")
  }
  
  fun <- "mrss.bcra4f3"
  par <- list(mdes=mdes, power=power, alpha=alpha, two.tail=two.tail,
              gm=gm, ncase=ncase, constrain=constrain,
              n=n, J=J, K=K, L=L, L0=L0, K0=K0, tol=tol, 
              rho2=rho2, rho3=rho3,
              P=P, g3=g3, R12=R12, R22=R22, R32=R32)
  mrss.out <- list(fun=fun,par=par,round.mrss=round(round.mrss, digits=3),
              integer.mrss=round(integer.mrss, digits=3))
  class(mrss.out) <- c("pars")
  print(round(round.mrss, digits=3))
  return(invisible(mrss.out))
}
  
# example
# mrss.bcra4f3(rho3=.15, rho2=.15, g3=1, P=.5, R32=.5, R22=.5, R12=.5, n=10, J=4, K=4)
# mrss.bcra4f3(rho3=.15, rho2=.15, g3=1, P=.5, R32=.5, R22=.5, R12=.5, n=10, J=4, L=14)
# mrss.bcra4f3(rho3=.15, rho2=.15, g3=1, P=.5, R32=.5, R22=.5, R12=.5, n=10, K=4, L=14)
# mrss.bcra4f3(rho3=.15, rho2=.15, g3=1, P=.5, R32=.5, R22=.5, R12=.5, J=4, K=4, L=14)


# require ("nloptr")
# nJKL0 = starting values, default: c(10,10,10,10)
# ncase = numer of optimal samples in the output
# gm = grid multiplier to search best integer solutions
# constrain = "mdes", "power", or "cost"
# optimizer = "auglag_cobyla","auglag_lbfgs", "auglag_mma", or "auglag_slsqp"
optimal.bcra4f3 <- function(cn, cJ, cK, cL, cost=NULL, n=NULL, J=NULL, K=NULL, L=NULL,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           nJKL0=c(10,10,10,10), ncase=10, gm=2, 
                           constrain="cost", optimizer="auglag_cobyla",
                           rho3, rho2, P=.50, g3=0,  R32=0, R22=0, R12=0){
  
  # equality constraints on cost
  eq.cost <- function(nJKL){
    n <- nJKL[1]
    J <- nJKL[2]
    K <- nJKL[3]
    L <- nJKL[4]
    cost0 <- cL*L + cK*K*L + cJ*J*K*L + cn*n*J*K*L - cost
    return(cost0)
  }
  
  # equality constraints on power
  eq.power <- function(nJKL){
    n <- nJKL[1]
    J <- nJKL[2]
    K <- nJKL[3]
    L <- nJKL[4]
    df <- (L*(K-2)-g3)[[1]]
    lamda <- mdes/sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                       rho2*(1-R22)/(P*(1-P)*J*K*L) +
                      (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
    power1 <- ifelse(two.tail==FALSE,
                     1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                     1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                       pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
    power0 <- power1-power
    return(power0)
  }
  
  # equality constraints on mdes
  eq.mdes <- function(nJKL){
    n <- nJKL[1]
    J <- nJKL[2]
    K <- nJKL[3]
    L <- nJKL[4]
    df <- (L*(K-2)-g3)[[1]]
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    mdes1 <- M*sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                   rho2*(1-R22)/(P*(1-P)*J*K*L) +
                  (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
    mdes0 <-mdes1-mdes
    return(mdes0)
  }

  # minimize cost given power or mdes
  min.cost <- function(nJKL){
    n <- nJKL[1]
    J <- nJKL[2]
    K <- nJKL[3]
    L <- nJKL[4]
    cost <- cL*L + cK*K*L + cJ*J*K*L + cn*n*J*K*L
    return(cost)
  }
  
  # minimize treatment variance given cost
  min.var <- function(nJKL){
    n <- nJKL[1]
    J <- nJKL[2]
    K <- nJKL[3]
    L <- nJKL[4]
    sigmaT2 <- rho3*(1-R32)/(P*(1-P)*K*L) +
               rho2*(1-R22)/(P*(1-P)*J*K*L) +
              (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n)
    return(sigmaT2)
  }
  
  fn.constr <- switch(constrain,
                      "power"= "eq.power",
                      "mdes"= "eq.mdes",
                      "cost"= "eq.cost"
  )
  
  fn.min <- switch(constrain,
                   "power"= "min.cost",
                   "mdes"= "min.cost",
                   "cost"= "min.var"
  )
  
  local.optim <- switch(optimizer,
                        "auglag_cobyla"="COBYLA",
                        "auglag_lbfgs"="LBFGS",
                        "auglag_mma"="MMA",
                        "auglag_slsqp"="SLSQP"
  )
  
  # constrain specified sample sizes, define limits
  n0  <- ifelse(!is.null(n), n, nJKL0[1])
  nLB <- ifelse(!is.null(n), n, 1)
  nUB <- ifelse(!is.null(n), n, Inf)
  J0  <- ifelse(!is.null(J), J, nJKL0[2])
  JLB <- ifelse(!is.null(J), J, 1)
  JUB <- ifelse(!is.null(J), J, Inf)
  K0  <- ifelse(!is.null(K), K, nJKL0[3])
  KLB <- ifelse(!is.null(K), K, 3)
  KUB <- ifelse(!is.null(K), K, Inf)
  L0  <- ifelse(!is.null(L), L, nJKL0[4])
  LLB <- ifelse(!is.null(L), L, 1+g3)
  LUB <- ifelse(!is.null(L), L, Inf)
  nJKL0  <- c(n0,J0,K0,L0)
  nJKLLB <- c(nLB,JLB,KLB,LLB)
  nJKLUB <- c(nUB,JUB,KUB,LUB)
  
  # given equality contraints find real numbers that minimizes the treatment variance or cost
  nlopt.nJKL <- nloptr::auglag(x0=nJKL0, fn=fn.min, heq=fn.constr, 
                               localsolver=local.optim, localtol=1e-8, 
                               lower=nJKLLB, upper=nJKLUB)
  if(nlopt.nJKL$convergence<0 | all(nlopt.nJKL$par==nJKL0)){warning("Solution may not be feasible. Change default settings.")}
  
  # round solution
  nJKL1 <- round(nlopt.nJKL$par)
  est.cost <-  min.cost(nJKL1)
  est.mdes <- eq.mdes(nJKL1) + mdes
  est.power <- eq.power(nJKL1) + power
  round.optim <- cbind(nJKL1[1], nJKL1[2], nJKL1[3], nJKL1[4], est.cost, est.mdes, est.power)
  colnames(round.optim) <- c("n", "J", "K", "L", "cost", "mdes", "power")

  # grid to approximate an integer solution
  if(gm<=0){stop("'gm' should be >= 1")}
  if(ncase<=1){stop("'ncase' should be >= 2")}
  gridnJKL <- expand.grid(c(nJKL1[1]-gm*4):(nJKL1[1]+gm*4),
                          c(nJKL1[2]-gm*3):(nJKL1[2]+gm*3),
                          c(nJKL1[3]-gm*2):(nJKL1[3]+gm*2),
                          c(nJKL1[4]-gm*1):(nJKL1[4]+gm*1))
  
  # ensure nonzero postive cells while considering cells associated with df
  gridnJKL <- as.matrix(gridnJKL[gridnJKL[,1]>0 &
                                 gridnJKL[,2]>0 &
                                 gridnJKL[,3]>2 &
                                 gridnJKL[,4]>g3, ])
  
  # difference from specified power (default, .80), mdes (default, .25), or cost
  integer.optim <- matrix(NA,nrow(gridnJKL),7)
  integer.optim[,1:4] <- gridnJKL
  for(i in 1:nrow(gridnJKL)){
    integer.optim[i,5] <- min.cost(gridnJKL[i,])
    integer.optim[i,6] <- eq.mdes(gridnJKL[i,])
    integer.optim[i,7] <- eq.power(gridnJKL[i,])
  }
  
  integer.optim <- round(integer.optim, digits=3)
  
  # output result
  if(constrain=="power"){
    idx <- order(abs(integer.optim[,7]), integer.optim[,5],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,6] <- integer.optim[,6] + mdes
    integer.optim[,7] <- integer.optim[,7] + power
    colnames(integer.optim) <- c("n", "J", "K", "L", "cost", "mdes", "power")
  }else if(constrain=="mdes"){
    idx <- order(abs(integer.optim[,6]), integer.optim[,5],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,6] <- integer.optim[,6] + mdes
    integer.optim[,7] <- integer.optim[,7] + power
    colnames(integer.optim) <- c("n", "J", "K", "L", "cost", "mdes", "power")
  }else if(constrain=="cost"){
    idx <- order(abs(integer.optim[,5]-cost), -integer.optim[,7],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,6] <- integer.optim[,6] + mdes
    integer.optim[,7] <- integer.optim[,7] + power
    colnames(integer.optim) <- c("n", "J", "K", "L", "cost", "mdes", "power")
  }
  
  fun <- "optimal.bcra4f3"
  par <- list(cn=cn, cJ=cJ, cK=cK, cL=cL, cost=cost, n=n, J=J, K=K, L=L,
              power=power, mdes=mdes, alpha=alpha, two.tail=two.tail,
              nJKL0=nJKL0, ncase=ncase, gm=gm, 
              constrain=constrain, optimizer=optimizer,
              rho2=rho2, rho3=rho3,
              P=P, g3=g3, R12=R12, R22=R22, R32=R32)
  nloptr <- list(pars=nlopt.nJKL$par,
                 obj.value=nlopt.nJKL$value,
                 iter=nlopt.nJKL$iter,
                 global.solver=nlopt.nJKL$global_solver,
                 local.solver=nlopt.nJKL$local_solver,
                 convergence=nlopt.nJKL$convergence,
                 message=nlopt.nJKL$message)
  optim.out <- list(fun=fun,
              par=par,
              nloptr=nloptr,
              round.optim=round(round.optim, digits=3),
              integer.optim=round(integer.optim, digits=3)
              )
  class(optim.out) <- c("pars")
  print(round(round.optim, digits=3))
  return(invisible(optim.out))
}

# examples

# optimal sample given total cost
# optimal.bcra4f3(cn=1, cJ=10, cK=100, cL=1000, cost=75600, constrain="cost", rho3=.15, rho2=.20)

# optimal sample given per unit costs and power 
# optimal.bcra4f3(cn=1, cJ=10, cK=100, cL=1000, power=.80, constrain="power", rho3=.15, rho2=.20, nJKL0=c(3,3,3,3))

# optimal sample given per unit costs and mdes
# optimal.bcra4f3(cn=1, cJ=10, cK=100, cL=1000, power=.80, constrain="mdes", rho3=.15, rho2=.20, nJKL0=c(3,3,3,3))

# conditional optimal sample (fixed n=10, and J=2) given total cost
# optimal.bcra4f3(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, cost=75600, constrain="cost", rho3=.15, rho2=.20)

# conditional optimal sample (fixed n=10, and J=2) given per unit costs and power 
# optimal.bcra4f3(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, power=.80, constrain="power", rho3=.15, rho2=.20, nJKL0=c(3,3,3,3))

# conditional optimal sample (fixed n=10, and J=2) given per unit costs and mdes
# optimal.bcra4f3(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, power=.80, constrain="mdes", rho3=.15, rho2=.20, nJKL0=c(3,3,3,3))

  
