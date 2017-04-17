mdes.bira2r1 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2,  omega2, P=.50, g2=0, R12=0, RT22=0,  
                        n, J, ...){
  df <- J-g2-1
  T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
  T2 <- abs(qt(power,df))
  M <- ifelse(power>=.5,T1+T2,T1-T2)
  mdes <- M*sqrt(rho2*omega2*(1-RT22)/J +
               (1-rho2)*(1-R12)/(P*(1-P)*J*n))
  LCI <- mdes*(1 - T1/M)
  UCI <- mdes*(1 + T1/M)
  MLU <- cbind(mdes, LCI, UCI)
  colnames(MLU) <- c("mdes", "95% LCL", "95% UCL")
  fun <- "mdes.bira2r1"
  par <- list(power=power, alpha=alpha, two.tail=two.tail, 
               rho2=rho2, omega2=omega2, P=P, g2=g2, R12=R12, RT22=RT22,
               n=n, J=J)
  mdes.out <- list(fun=fun,par=par,df=df,M=M,mdes=round(MLU,digits=3))
  class(mdes.out) <- c("pars")
  print(round(MLU, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.bira2r1(rho2=.35, omega2=.10, n=83, J=480)


power.bira2r1 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2,  omega2, g2=0, P=.50, R12=0, RT22=0,
                         n, J, ...){
  df <- J-g2-1
  lamda <- mdes/sqrt(rho2*omega2*(1-RT22)/J +
                    (1-rho2)*(1-R12)/(P*(1-P)*J*n))
  power <- ifelse(two.tail==FALSE,
                  1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                  1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                    pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
  
  fun <- "power.bira2r1"
  par <- list(mdes=mdes, alpha=alpha, two.tail=two.tail, 
              rho2=rho2, omega2=omega2, P=P, g2=g2, R12=R12, RT22=RT22,
              n=n, J=J)
  power.out <- list(fun=fun,par=par,df=df,lamda=lamda,power=round(power,digits=3))
  class(power.out) <- c("pars")
  print(paste("power = ", round(power, digits=3)))
  return(invisible(power.out))
}
# example
# power.bira2r1(rho2=.35, omega2=.10, n=83, J=480)

# gm = grid multiplier: default 2
# ncase = numer of fixed samples in the output
# constrain = "power", or "mdes" :  default "power"
mrss.bira2r1 <- function(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
                        gm=2, ncase=10, constrain="power",
                        n=NULL, J=NULL, J0=10, tol=.10,
                        rho2, omega2, g2=0, P=.50, R12=0, RT22=0){
  
  # sample size for one level allowed to be null, the rest must be specified
  check.NULL <- c(is.null(n),is.null(J))
  if(length(check.NULL[check.NULL==TRUE])>1){
    stop("Specify one of the n, J")
  }
  
  # equality constraints on power
  eq.power <- function(nJ){
    n <- nJ[1]
    J <- nJ[2]
    df <- (J-g2-1)[[1]]
    lamda <- mdes/sqrt(rho2*omega2*(1-RT22)/J +
                      (1-rho2)*(1-R12)/(P*(1-P)*J*n))
    power1 <- ifelse(two.tail==FALSE,
                    1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                    1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                      pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
    power0 <- power1-power
    return(power0)
  }
  
  # equality constraints on mdes
  eq.mdes <- function(nJ){
    n <- nJ[1]
    J <- nJ[2]
    df <- (J-g2-1)[[1]]
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    mdes1 <- M*sqrt(rho2*omega2*(1-RT22)/J +
                  (1-rho2)*(1-R12)/(P*(1-P)*J*n))
    mdes0 <-mdes1-mdes
    return(mdes0)
  }
  
  # find minimum required sample size per higher unit
  if(is.null(J)){
    i <- 0
    conv <- FALSE
    while(i<=100 & conv==FALSE){
      df <- J0-g2-1
      if(df<= 0 | is.infinite(df)){break}
      T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
      T2 <- abs(qt(power,df))
      M <- ifelse(power>=.5,T1+T2,T1-T2)
      J1 <- (M/mdes)^2 * (rho2*omega2*(1-RT22) +
                            (1-rho2)*(1-R12)/(P*(1-P)*n))
      if(abs(J1-J0)<tol){conv <- TRUE}
      J0 <- (J1+J0)/2
      i <- i+1
    }
    J <- ifelse(df>0,round(J0),NA)
  }else if(is.null(n)){
    df <- J-g2-1
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    n <- (1-rho2)*(1-R12) /
        (P*(1-P)*(J*(mdes/M)^2 - rho2*omega2*(1-RT22)))
    n <- round(n)
  }
  
  nJ1 <- c(n,J)
  if(any(nJ1<=0)|any(is.na(nJ1))){stop("Solution not feasible due to nonpositive sample size estimation.")}
  
  # round solution
  est.mdes <- eq.mdes(nJ1) + mdes
  est.power <- eq.power(nJ1) + power
  round.mrss <- cbind(nJ1[1], nJ1[2], est.mdes, est.power)
  colnames(round.mrss) <- c("n", "J", "mdes", "power")
 
  # grid to approximate an integer solution
  if(gm<=0){stop("'gm' should be >= 1")}
  if(ncase<=1){stop("'ncase' should be >= 2")}
  gridnJ <- expand.grid(c(nJ1[1]-gm*2):(nJ1[1]+gm*2),
                        c(nJ1[2]-gm*1):(nJ1[2]+gm*1))
  
  # ensure nonzero postive cells while considering cells associated with df
  gridnJ <- as.matrix(gridnJ[gridnJ[,1]>0 &
                             gridnJ[,2]>g2+1, ])
  
  # difference from specified power (default, .80), mdes (default, .25), or cost
  integer.mrss <- matrix(NA,nrow(gridnJ),4)
  integer.mrss[,1:2] <- gridnJ
  for(i in 1:nrow(gridnJ)){
    integer.mrss[i,3] <- eq.mdes(gridnJ[i,])
    integer.mrss[i,4] <- eq.power(gridnJ[i,])
  }
  
  # output result
  if(constrain=="power"){
    idx <- order(abs(integer.mrss[,4]),
                 decreasing=FALSE)[1:ncase]
    integer.mrss <- integer.mrss[idx,]
    integer.mrss[,3] <- integer.mrss[,3] + mdes
    integer.mrss[,4] <- integer.mrss[,4] + power
    colnames(integer.mrss) <- c("n", "J","mdes", "power")
  }else if(constrain=="mdes"){
    idx <- order(abs(integer.mrss[,3]),
                 decreasing=FALSE)[1:ncase]
    integer.mrss <- integer.mrss[idx,]
    integer.mrss[,3] <- integer.mrss[,3] + mdes
    integer.mrss[,4] <- integer.mrss[,4] + power
    colnames(integer.mrss) <- c("n", "J", "mdes", "power")
  }
  
  fun <- "mrss.bira2r1"
  par <- list(mdes=mdes, power=power, alpha=alpha, two.tail=two.tail,
              gm=gm, ncase=ncase, constrain=constrain,
              n=n, J=J, J0=J0, tol=tol,
              rho2=rho2, omega2=omega2, g2=g2, P=P, R12=R12, RT22=RT22)
  mrss.out <- list(fun=fun,par=par,round.mrss=round(round.mrss, digits=3),
              integer.mrss=round(integer.mrss, digits=3))
  class(mrss.out) <- c("pars")
  print(round(round.mrss, digits=3))
  return(invisible(mrss.out))
}
  
# example
# mrss.bira2r1(rho2=.35, omega2=.10, n=83)
# mrss.bira2r1(rho2=.35, omega2=.10, J=10) 

# require ("nloptr")
# nJ0 = starting values, default: c(10,10)
# ncase = numer of optimal samples in the output
# gm = grid multiplier to search best integer solutions
# constrain = "mdes", "power", or "cost"
# optimizer = "auglag_cobyla","auglag_lbfgs", "auglag_mma", or "auglag_slsqp"
optimal.bira2r1 <- function(cn, cJ, cost=NULL, n=NULL, J=NULL,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           nJ0=c(10,10), ncase=10, gm=2, 
                           constrain="cost", optimizer="auglag_cobyla",
                           rho2, omega2, P=.50, g2=0, R12=0, RT22=0){
  
  # equality constraints on cost
  eq.cost <- function(nJ){
    n <- nJ[1]
    J <- nJ[2]
    cost0 <- cJ*J + cn*n*J - cost
    return(cost0)
  }
  
  # equality constraints on power
  eq.power <- function(nJ){
    n <- nJ[1]
    J <- nJ[2]
    df <- (J-g2-1)[[1]]
    lamda <- mdes/sqrt(rho2*omega2*(1-RT22)/J +
                      (1-rho2)*(1-R12)/(P*(1-P)*J*n))
    power1 <- ifelse(two.tail==FALSE,
                     1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                     1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                       pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
    power0 <- power1-power
    return(power0)
  }
  
  # equality constraints on mdes
  eq.mdes <- function(nJ){
    n <- nJ[1]
    J <- nJ[2]
    df <- (J-g2-1)[[1]]
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    mdes1 <- M*sqrt(rho2*omega2*(1-RT22)/J +
                  (1-rho2)*(1-R12)/(P*(1-P)*J*n))
    mdes0 <-mdes1-mdes
    return(mdes0)
  }

  # minimize cost given power or mdes
  min.cost <- function(nJ){
    n <- nJ[1]
    J <- nJ[2]
    cost <- cJ*J + cn*n*J
    return(cost)
  }
  
  # minimize treatment variance given cost
  min.var <- function(nJ){
    n <- nJ[1]
    J <- nJ[2]
    sigmaT2 <- rho2*omega2*(1-RT22)/J +
               (1-rho2)*(1-R12)/(P*(1-P)*J*n)
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
  n0  <- ifelse(!is.null(n), n, nJ0[1])
  nLB <- ifelse(!is.null(n), n, 1)
  nUB <- ifelse(!is.null(n), n, Inf)
  J0  <- ifelse(!is.null(J), J, nJ0[2])
  JLB <- ifelse(!is.null(J), J, g2+2)
  JUB <- ifelse(!is.null(J), J, Inf)
  nJ0  <- c(n0,J0)
  nJLB <- c(nLB,JLB)
  nJUB <- c(nUB,JUB)
  
  # given equality contraints find real numbers that minimizes the treatment variance or cost
  nlopt.nJ <- nloptr::auglag(x0=nJ0, fn=fn.min, heq=fn.constr, 
                               localsolver=local.optim, localtol=1e-8, 
                               lower=nJLB, upper=nJUB)
  if(nlopt.nJ$convergence<0 | all(nlopt.nJ$par==nJ0)){warning("Solution may not be feasible. Change default settings.")}
  
  # round solution
  nJ1 <- round(nlopt.nJ$par)
  est.cost <-  min.cost(nJ1)
  est.mdes <- eq.mdes(nJ1) + mdes
  est.power <- eq.power(nJ1) + power
  round.optim <- cbind(nJ1[1], nJ1[2], est.cost, est.mdes, est.power)
  colnames(round.optim) <- c("n", "J","cost", "mdes", "power")
  
  # grid to approximate an integer solution
  if(gm<=0){stop("'gm' should be >= 1")}
  if(ncase<=1){stop("'ncase' should be >= 2")}
  gridnJ <- expand.grid(c(nJ1[1]-gm*2):(nJ1[1]+gm*2),
                        c(nJ1[2]-gm*1):(nJ1[2]+gm*1))
  
  # ensure nonzero postive cells while considering cells associated with df
  gridnJ <- as.matrix(gridnJ[gridnJ[,1]>0 &
                             gridnJ[,2]>g2+1, ])
  
  # difference from specified power (default, .80), mdes (default, .25), or cost
  integer.optim <- matrix(NA,nrow(gridnJ),5)
  integer.optim[,1:2] <- gridnJ
  for(i in 1:nrow(gridnJ)){
    integer.optim[i,3] <- min.cost(gridnJ[i,])
    integer.optim[i,4] <- eq.mdes(gridnJ[i,])
    integer.optim[i,5] <- eq.power(gridnJ[i,])
  }
  
  integer.optim <- round(integer.optim, digits=3)
  
  # output result
  if(constrain=="power"){
    idx <- order(abs(integer.optim[,5]), integer.optim[,3],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,4] <- integer.optim[,4] + mdes
    integer.optim[,5] <- integer.optim[,5] + power
    colnames(integer.optim) <- c("n", "J", "cost", "mdes", "power")
  }else if(constrain=="mdes"){
    idx <- order(abs(integer.optim[,4]), integer.optim[,3],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,4] <- integer.optim[,4] + mdes
    integer.optim[,5] <- integer.optim[,5] + power
    colnames(integer.optim) <- c("n", "J", "cost", "mdes", "power")
  }else if(constrain=="cost"){
    idx <- order(abs(integer.optim[,3]-cost), -integer.optim[,5],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,4] <- integer.optim[,4] + mdes
    integer.optim[,5] <- integer.optim[,5] + power
    colnames(integer.optim) <- c("n", "J", "cost", "mdes", "power")
  }
  
  fun <- "optimal.bira2r1"
  par <- list(cn=cn, cJ=cJ, cost=cost, n=n, J=J,
              power=power, mdes=mdes, alpha=alpha, two.tail=two.tail,
              nJ0=nJ0, ncase=ncase, gm=gm, 
              constrain=constrain, optimizer=optimizer,
              rho2=rho2, omega2=omega2, P=P, g2=g2, R12=R12, RT22=RT22)
  nloptr <- list(pars=nlopt.nJ$par,
                 obj.value=nlopt.nJ$value,
                 iter=nlopt.nJ$iter,
                 global.solver=nlopt.nJ$global_solver,
                 local.solver=nlopt.nJ$local_solver,
                 convergence=nlopt.nJ$convergence,
                 message=nlopt.nJ$message)
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
# optimal.bira2r1(cn=1, cJ=10, cost=560, constrain="cost", rho2=.20, omega2=.50)

# optimal sample given per unit costs and power 
# optimal.bira2r1(cn=1, cJ=10, power=.80, constrain="power", rho2=.20, omega2=.50)

# optimal sample given per unit costs and mdes
# optimal.bira2r1(cn=1, cJ=10, power=.80, constrain="mdes", rho2=.20, omega2=.50)

# conditional optimal sample (fixed n=10) given total cost
# optimal.bira2r1(cn=1, cJ=10, n=10, cost=560, constrain="cost", rho2=.20, omega2=.50)

# conditional optimal sample (fixed n=10) given per unit costs and power 
# optimal.bira2r1(cn=1, cJ=10, n=10, constrain="power", rho2=.20, omega2=.50)

# conditional optimal sample (fixed n=10) given per unit costs and mdes
# optimal.bira2r1(cn=1, cJ=10, n=10, constrain="mdes", rho2=.20, omega2=.50)

  
