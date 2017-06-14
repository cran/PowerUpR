mdes.ira1r1 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        P=.50, g1=0, R12=0, n, ...){
  df <- n-g1-2
  T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
  T2 <- abs(qt(power,df))
  M <- ifelse(power>=.5,T1+T2,T1-T2)
  mdes <- M*sqrt((1-R12)/(P*(1-P)*n))
  LCI <- mdes*(1 - T1/M)
  UCI <- mdes*(1 + T1/M)
  MLU <- cbind(mdes, LCI, UCI)
  colnames(MLU) <- c("mdes", "95% LCL", "95% UCL")
  fun <- "mdes.ira1r1"
  par <- list(power=power, alpha=alpha, two.tail=two.tail, 
              P=P, g1=g1, R12=R12, n=n)
  mdes.out <- list(fun=fun,par=par,df=df,M=M,mdes=round(MLU,digits=3))
  class(mdes.out) <- c("pars")
  print(round(MLU, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.ira1r1(n=400)

power.ira1r1 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                        P=.50, g1=0, R12=0, n, ...){
  df <- n-g1-2
  lamda <- mdes/sqrt((1-R12)/(P*(1-P)*n))
  power <- ifelse(two.tail==FALSE,
                  1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                  1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                    pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
  fun <- "power.ira1r1"
  par <- list(mdes=mdes, alpha=alpha, two.tail=two.tail, 
              P=P, g1=g1, R12=R12, n=n)
  power.out <- list(fun=fun,par=par,df=df,lamda=lamda,power=round(power,digits=3))
  class(power.out) <- c("pars")
  print(paste("power = ", round(power, digits=3)))
  return(invisible(power.out))
}
# example
# power.ira1r1(n=400)

# gm = grid multiplier: default 2
# ncase = numer of fixed samples in the output
# constrain = "power", or "mdes" :  default "power"
mrss.ira1r1 <- function(mdes=.25, power=.80, alpha=.05, two.tail=TRUE,
                        gm=10, ncase=10, constrain="power",
                        n=NULL, n0=10, tol=.10,
                        P=.50, g1=0, R12=0){
  
  # equality constraints on power
  eq.power <- function(N){
    n <- N[1]
    df <- (n-g1-2)[[1]]
    lamda <- mdes/sqrt((1-R12)/(P*(1-P)*n))
    power1 <- ifelse(two.tail==FALSE,
                    1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                    1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                      pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
    power0 <- power1-power
    return(power0)
  }
  
  # equality constraints on mdes
  eq.mdes <- function(N){
    n <- N[1]
    df <- (n-g1-2)[[1]]
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    mdes1 <- M*sqrt((1-R12)/(P*(1-P)*n))
    mdes0 <-mdes1-mdes
    return(mdes0)
  }
  
  # find minimum required sample size 
  if(is.null(n)){
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
    n <- ifelse(df>0,round(n0),NA)
  }
  
  N1 <- c(n)
  if(any(N1<=0)|any(is.na(N1))){stop("Solution not feasible due to nonpositive sample size estimation.")}
  
  # round solution
  est.mdes <- eq.mdes(N1) + mdes
  est.power <- eq.power(N1) + power
  round.mrss <- cbind(N1[1], est.mdes, est.power)
  colnames(round.mrss) <- c("n", "mdes", "power")
 
  # grid to approximate an integer solution 
  if(gm<=0){stop("'gm' should be >= 1")}
  if(ncase<=1){stop("'ncase' should be >= 2")}
  gridN <- expand.grid(c(N1[1]-gm):(N1[1]+gm))
  
  # ensure nonzero postive cells while considering cells associated with df
  gridN <- as.matrix(gridN[gridN[,1]>g1+2, ])
  
  # difference from specified power (default, .80), mdes (default, .25), or cost
  integer.mrss <- matrix(NA,nrow(gridN),3)
  integer.mrss[,1] <- gridN
  for(i in 1:nrow(gridN)){
    integer.mrss[i,2] <- eq.mdes(gridN[i,])
    integer.mrss[i,3] <- eq.power(gridN[i,])
  }
  
  # output result
  if(constrain=="power"){
    idx <- order(abs(integer.mrss[,3]),
                 decreasing=FALSE)[1:ncase]
    integer.mrss <- integer.mrss[idx,]
    integer.mrss[,2] <- integer.mrss[,2] + mdes
    integer.mrss[,3] <- integer.mrss[,3] + power
    colnames(integer.mrss) <- c("n", "mdes", "power")
  }else if(constrain=="mdes"){
    idx <- order(abs(integer.mrss[,2]),
                 decreasing=FALSE)[1:ncase]
    integer.mrss <- integer.mrss[idx,]
    integer.mrss[,2] <- integer.mrss[,2] + mdes
    integer.mrss[,3] <- integer.mrss[,3] + power
    colnames(integer.mrss) <- c("n", "mdes", "power")
  }
  fun <- "mrss.ira1r1"
  par <- list(mdes=mdes, power=power, alpha=alpha, two.tail=two.tail,
              gm=gm, ncase=ncase, constrain=constrain,
              n=n, n0=n0, tol=tol, 
              P=P, g1=g1, R12=R12)
  mrss.out <- list(fun=fun,par=par,round.mrss=round(round.mrss, digits=3),
              integer.mrss=round(integer.mrss, digits=3))
  class(mrss.out) <- c("pars")
  print(round(round.mrss, digits=3))
  return(invisible(mrss.out))
}
  
# example
# mrss.ira1r1()

# require ("nloptr")
# N0 = starting values, default: c(10)
# ncase = numer of optimal samples in the output
# gm = grid multiplier to search best integer solutions
# constrain = "mdes", "power", or "cost"
# optimizer = "auglag_cobyla","auglag_lbfgs", "auglag_mma", or "auglag_slsqp"
optimal.ira1r1 <- function(cn, cost=NULL, n=NULL,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           N0=c(10), ncase=10, gm=10, 
                           constrain="cost", optimizer="auglag_cobyla",
                           P=.50, g1=0, R12=0){
  
  # equality constraints on cost
  eq.cost <- function(N){
    n <- N[1]
    cost0 <- cn*n - cost
    return(cost0)
  }
  
  # equality constraints on power
  eq.power <- function(N){
    n <- N[1]
    df <- (n-g1-2)[[1]]
    lamda <- mdes/sqrt((1-R12)/(P*(1-P)*n))
    power1 <- ifelse(two.tail==FALSE,
                     1-pt(qt(alpha,df,lower.tail=FALSE),df,lamda),
                     1-pt(qt(alpha/2,df,lower.tail=FALSE),df,lamda)+
                       pt(-qt(alpha/2,df,lower.tail=FALSE),df,lamda))
    power0 <- power1-power
    return(power0)
  }
  
  # equality constraints on mdes
  eq.mdes <- function(N){
    n <- N[1]
    df <- (n-g1-2)[[1]]
    T1 <- ifelse(two.tail==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    mdes1 <- M*sqrt((1-R12)/(P*(1-P)*n))
    mdes0 <-mdes1-mdes
    return(mdes0)
  }

  # minimize cost given power or mdes
  min.cost <- function(N){
    n <- N[1]
    cost <- cn*n
    return(cost)
  }
  
  # minimize treatment variance given cost
  min.var <- function(N){
    n <- N[1]
    sigmaT2 <- (1-R12)/(P*(1-P)*n)
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
  n0  <- ifelse(!is.null(n), n, N0[1])
  nLB <- ifelse(!is.null(n), n, g1+2+1)
  nUB <- ifelse(!is.null(n), n, Inf)
  N0  <- c(n0)
  NLB <- c(nLB)
  NUB <- c(nUB)
  
  # given equality contraints find real numbers that minimizes the treatment variance or cost
  nlopt.N <- nloptr::auglag(x0=N0, fn=fn.min, heq=fn.constr, 
                               localsolver=local.optim, localtol=1e-8, 
                               lower=NLB, upper=NUB)
  if(nlopt.N$convergence<0 | all(nlopt.N$par==N0)){warning("Solution may not be feasible. Change default settings.")}
  
  # round solution
  N1 <- round(nlopt.N$par)
  est.cost <-  min.cost(N1)
  est.mdes <- eq.mdes(N1) + mdes
  est.power <- eq.power(N1) + power
  round.optim <- cbind(N1[1], est.cost, est.mdes, est.power)
  colnames(round.optim) <- c("n", "cost", "mdes", "power")

  # estimated solution
  est.cost <-  min.cost(N1)
  est.mdes <- eq.mdes(N1) + mdes
  est.power <- eq.power(N1) + power
  est.optim <- cbind(N1[1], est.cost, est.mdes, est.power)
  colnames(est.optim) <- c("n", "cost", "mdes", "power")
  
  # grid to approximate an integer solution
  if(gm<=0){stop("'gm' should be >= 1")}
  if(ncase<=1){stop("'ncase' should be >= 2")}
  gridN <- expand.grid(c(N1[1]-gm):(N1[1]+gm))
  
  # ensure nonzero postive cells while considering cells associated with df
  gridN <- as.matrix(gridN[gridN[,1]>g1+2, ])
  
  # difference from specified power (default, .80), mdes (default, .25), or cost
  integer.optim <- matrix(NA,nrow(gridN),4)
  integer.optim[,1] <- gridN
  for(i in 1:nrow(gridN)){
    integer.optim[i,2] <- min.cost(gridN[i,])
    integer.optim[i,3] <- eq.mdes(gridN[i,])
    integer.optim[i,4] <- eq.power(gridN[i,])
  }
  
  integer.optim <- round(integer.optim, digits=3)
  
  # output result
  if(constrain=="power"){
    idx <- order(abs(integer.optim[,4]), integer.optim[,2],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,3] <- integer.optim[,3] + mdes
    integer.optim[,4] <- integer.optim[,4] + power
    colnames(integer.optim) <- c("n", "cost", "mdes", "power")
  }else if(constrain=="mdes"){
    idx <- order(abs(integer.optim[,3]), integer.optim[,2],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,3] <- integer.optim[,3] + mdes
    integer.optim[,4] <- integer.optim[,4] + power
    colnames(integer.optim) <- c("n", "cost", "mdes", "power")
  }else if(constrain=="cost"){
    idx <- order(abs(integer.optim[,2]-cost), -integer.optim[,4],
                 decreasing=FALSE)[1:ncase]
    integer.optim <- integer.optim[idx,]
    integer.optim[,3] <- integer.optim[,3] + mdes
    integer.optim[,4] <- integer.optim[,4] + power
    colnames(integer.optim) <- c("n", "cost", "mdes", "power")
  }
  fun <- "optimal.ira1r1"
  par <- list(cn=cn, cost=cost, n=n,
              power=power, mdes=mdes, alpha=alpha, two.tail=two.tail,
              N0=N0, ncase=ncase, gm=gm, 
              constrain=constrain, optimizer=optimizer,
              P=P, g1=g1, R12=R12)
  nloptr <- list(pars=nlopt.N$par,
                 obj.value=nlopt.N$value,
                 iter=nlopt.N$iter,
                 global.solver=nlopt.N$global_solver,
                 local.solver=nlopt.N$local_solver,
                 convergence=nlopt.N$convergence,
                 message=nlopt.N$message)
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
# optimal.ira1r1(cn=1,  cost=560, constrain="cost")

# optimal sample given per unit costs and power 
# optimal.ira1r1(cn=1,  constrain="power")

# optimal sample given per unit costs and mdes
# optimal.ira1r1(cn=1,  constrain="mdes")


  
