globalVariables("ncp")

# normal (Schochet, 2008, p.14)
# P: proportion of cases in treatment
.RTZ.n <- function(P){
  RTZ <- dnorm(qnorm(1-P)) / sqrt(P*(1-P))
  return(RTZ)
}

# uniform (Schochet, 2008, p.14)
# P: proportion of cases in treatment
.RTZ.u <- function(P){
  RTZ <- sqrt(3*P*(1-P))
  return(RTZ)
}

# truncated normal (Schochet, 2008, p.14)
# P: proportion of cases in treatment
# k1: left truncation point (in standard deviation units from full normal distribution mean)
# k2: right truncation point (in standard deviation units from full normal distribution mean)
.RTZ.tn <- function(k1, k2, P){
  c <- qnorm(P*pnorm(k1) + (1-P)*pnorm(k2) )
  sigmaZ2 <- 1 - ((k2*dnorm(k2) - k1*dnorm(k1)) / (pnorm(k2) - pnorm(k1))) -
    ((dnorm(k2) - dnorm(k1))/(pnorm(k2) - pnorm(k1)))^2
  RTZ <- (P/sqrt(sigmaZ2*P*(1-P))) * ((dnorm(k2) - dnorm(k1)) / (pnorm(k2) - pnorm(k1)) -
                                        (dnorm(k2) - dnorm(c)) / (pnorm(k2) - pnorm(c)))
  return(RTZ)
}

# design effect
.D.fun <- function(P, k1, k2, dist.Z, RTZ){
  .DE <- function(RTZ){
    D <- 1/(1-RTZ^2)
    return(D)
  }
  if(is.null(RTZ)){
    if(dist.Z=="normal"){
      D <- .DE(.RTZ.tn(k1,k2,P))
    }else if(dist.Z=="uniform"){
      if(k1!=-6 | k2!=6){warning("k1 and/or k2 will be ignored.")}
      D <- .DE(.RTZ.u(P))
    }
  }else if(!is.null(RTZ)){
    D <- .DE(RTZ)
    warning("Make sure RTZ is consistent with P as both are related! See Schochet (2008, p.14).
            With RTZ=0 results are identical to random assignment case.")
  }
  return(D)
  }

# minimum detectable effect size
.mdes.fun <- function(power, alpha, SSE, df, two.tail, envir = parent.frame()){
  T1 <- ifelse(two.tail == TRUE, abs(qt(alpha/2, df)), abs(qt(alpha, df)))
  T2 <- abs(qt(power, df))
  M <- envir$ncp <- ifelse(power >= 0.5, T1 + T2, T1 - T2)
  mdes <- M * SSE
  LCI <- mdes * (1 - T1/M)
  UCI <- mdes * (1 + T1/M)
  MLU <- cbind(mdes, LCI, UCI)
  colnames(MLU) <- c("mdes", paste0(100 * (1 - alpha), "% LCL"), paste0(100 * (1 - alpha), "% LCL"))
  return(MLU)
}

# statistical power
.power.fun <- function(mdes, alpha, SSE, df, two.tail, envir = parent.frame()){
  lamda <- envir$ncp <- mdes/SSE
  power <- ifelse(two.tail == FALSE,
                  1 - pt(qt(alpha, df, lower.tail = FALSE), df, lamda),
                  1 - pt(qt(alpha/2, df, lower.tail = FALSE), df, lamda) +
                    pt(-qt(alpha/2, df, lower.tail = FALSE), df, lamda))
  return(power)
}

# constrained optimal sample allocation
.cosa.fun <- function(cn = 0, cJ = 0, cK = 0, cL = 0, cost = NULL,
                      n = NULL, J = NULL, K = NULL, L = NULL, P = NULL,
                      nJ0 = c(10,10), nJK0 = c(10,10,10), nJKL0 = c(10,10,10,10), P0 = .50,
                      ncase = 10, gm = 2, constrain = "power", optimizer = "auglag_cobyla",
                      rho2, rho3, rho4, omega2, omega3, omega4,
                      Q=NULL, g1=0, g2=0, g3=0, g4=0,
                      R12=0, R22=0, R32=0, R42=0, RT22=0, RT32=0, RT42=0, R2T2=0, R3T2=0,
                      ...){
  fun <- get("fun", parent.frame())
  df <- get("df", parent.frame())
  SSE <- get("SSE", parent.frame())
  two.tail <- get("two.tail", parent.frame())
  alpha <- get("alpha", parent.frame())
  power <- get("power", parent.frame())
  mdes <- get("mdes", parent.frame())
  LB <- get("LB", parent.frame())
  nlevels <- length(LB)

  # check for valid unequal costs, which applies to levels at or below randomized level
  rlevel <- as.numeric(substr(fun, nchar(fun), nchar(fun)))
  if(nlevels==2 & rlevel==1 & length(cJ)>2){
    stop(simpleError(paste(sQuote(quote(cJ)), "cannot have a length greater than 1")))
  }
  if(nlevels==3){
    if(rlevel==1){
      if(length(cJ)==2){
        stop(simpleError(paste(sQuote(quote(cJ)), "cannot have a length greater than 1")))
      }
      if(length(cK)==2){
        stop(simpleError(paste(sQuote(quote(cK)), "cannot have a length greater than 1")))
      }
    }else if(rlevel==2 & length(cK)==2){
      stop(simpleError(paste(sQuote(quote(cK)), "cannot have a length greater than 1")))
    }
  }
  if(nlevels==4){
    if(rlevel==1){
      if(length(cJ)==2){
        stop(simpleError(paste(sQuote(quote(cJ)), "cannot have a length greater than 1")))
      }else if(length(cK)==2){
        stop(simpleError(paste(sQuote(quote(cK)), "cannot have a length greater than 1")))
      }else if(length(cL)==2){
        stop(simpleError(paste(sQuote(quote(cL)), "cannot have a length greater than 1")))
      }
    }else if(rlevel==3 & length(cL)==2){
      stop(simpleError(paste(sQuote(quote(cL)), "cannot have a length greater than 1")))
    }
  }

  # check for valid cost information
  if(length(cn)==1){
    cn <- c(cn,cn)
  }else if(length(cn)==2){
    warning(simpleWarning(paste("The first value in", sQuote(quote(cn)), "is the cost associated with treatment units" )))
  }else if(length(cn)>2){
    stop(simpleError(paste(sQuote(quote(cn)), "cannot have a length greater than 2")))
  }

  if(length(cJ)==1){
    cJ <- c(cJ,cJ)
  }else if(length(cJ)==2){
    warning(simpleWarning(paste("The first value in", sQuote(quote(cJ)), "is the cost associated with treatment units" )))
  }else if(length(cJ)>2){
    stop(simpleError(paste(sQuote(quote(cJ)), "cannot have a length greater than 2")))
  }

  if(length(cK)==1){
    cK <- c(cK,cK)
  }else if(length(cK)==2){
    warning(simpleWarning(paste("The first value in", sQuote(quote(cK)), "is the cost associated with treatment units" )))
  }else if(length(cK)>2){
    stop(simpleError(paste(sQuote(quote(cK)), "cannot have a length greater than 2")))
  }

  if(length(cL)==1){
    cL <- c(cL,cL)
  }else if(length(cL)==2){
    warning(simpleWarning(paste("The first value in", sQuote(quote(cL)), "is the cost associated with treatment units" )))
  }else if(length(cL)>2){
    stop(simpleError(paste(sQuote(quote(cL)), "cannot have a length greater than 2")))
  }

  if(!is.null(cost)){
    if(length(cost)>1){
      stop(simpleError(paste(sQuote(quote(cost)), "cannot have a length greater than 1")))
    }else if(cost<0){
      stop(simpleError(paste(sQuote(quote(cost)), "cannot have a negative value")))
    }
  }

  # check for valid P in RDD designs (not under researcher control)
  fun.parsed <- scan(text = fun, what = "character", sep=".", quiet = TRUE)
  if(fun.parsed[2] %in% c("crd2r2", "crd3r3", "bird2r1", "bird2f1", "bcrd3r2", "bcrd3f2")){
    if(length(P)>1 | !is.null(P)){
      stop(simpleError(paste(sQuote(quote(P)), "cannot be NULL or have a length greater than 1")))
    }
  }

  # equality constraints on cost
  .eq.cost <- function(ssP){
    n <- ssP[1]
    J <- ssP[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ssP[3]
    }
    if(nlevels == 4){
      L <- ssP[4]
    }
    P <- ssP[nlevels+1]
    cost0 <- switch(as.character(nlevels),
                    '2' = {P*(cJ[1]*J + cn[1]*n*J) + (1-P)*(cJ[2]*J + cn[2]*n*J) - cost},
                    '3' = {P*(cK[1]*K + cJ[1]*J*K + cn[1]*n*J*K) + (1-P)*(cK[2]*K + cJ[2]*J*K + cn[2]*n*J*K) - cost},
                    '4' = {P*(cL[1]*L + cK[1]*K*L + cJ[1]*J*K*L + cn[1]*n*J*K*L) + (1-P)*(cL[2]*L + cK[2]*K*L + cJ[2]*J*K*L + cn[2]*n*J*K*L) - cost})
    return(cost0)
  }

  # equality constraints on power
  .eq.power <- function(ssP){
    n <- ssP[1]
    J <- ssP[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ssP[3]
    }
    if(nlevels == 4){
      L <- ssP[4]
    }
    P <- ssP[nlevels+1]
    power1 <- .power.fun(mdes, alpha, eval(SSE), eval(df), two.tail)
    power0 <- power1-power
    return(power0)
  }

  # equality constraints on mdes
  .eq.mdes <- function(ssP){
    n <- ssP[1]
    J <- ssP[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ssP[3]
    }
    if(nlevels == 4){
      L <- ssP[4]
    }
    P <- ssP[nlevels+1]
    mdes1 <- .mdes.fun(power, alpha, eval(SSE), eval(df), two.tail)[1]
    mdes0 <-mdes1-mdes
    return(mdes0)
  }

  # minimize cost given power or mdes
  .min.cost <- function(ssP){
    n <- ssP[1]
    J <- ssP[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ssP[3]
    }
    if(nlevels == 4){
      L <- ssP[4]
    }
    P <- ssP[nlevels+1]
    cost <- switch(as.character(nlevels),
                   '2' = {P*(cJ[1]*J + cn[1]*n*J) + (1-P)*(cJ[2]*J + cn[2]*n*J)},
                   '3' = {P*(cK[1]*K + cJ[1]*J*K + cn[1]*n*J*K) + (1-P)*(cK[2]*K + cJ[2]*J*K + cn[2]*n*J*K)},
                   '4' = {P*(cL[1]*L + cK[1]*K*L + cJ[1]*J*K*L + cn[1]*n*J*K*L) + (1-P)*(cL[2]*L + cK[2]*K*L + cJ[2]*J*K*L + cn[2]*n*J*K*L)})
    return(cost)
  }

  # minimize treatment variance given cost
  .min.var <- function(ssP){
    n <- ssP[1]
    J <- ssP[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ssP[3]
    }
    if(nlevels == 4){
      L <- ssP[4]
    }
    P <- ssP[nlevels+1]
    sigmaT2 <- eval(SSE)^2
    return(sigmaT2)
  }

  fn.constr <- switch(constrain,
                      "power" = .eq.power,
                      "mdes" = .eq.mdes,
                      "cost" = .eq.cost)
  fn.min <- switch(constrain,
                   "power" = .min.cost,
                   "mdes" = .min.cost,
                   "cost" = .min.var)
  local.optim <- switch(optimizer,
                        "auglag_cobyla" = "COBYLA",
                        "auglag_lbfgs" = "LBFGS",
                        "auglag_mma" = "MMA",
                        "auglag_slsqp" = "SLSQP")

  # constraints on one or more of n, J, K, L and P
  P0  <- ifelse(!is.null(P), mean(P), P0)
  PLB <- ifelse(!is.null(P), min(P), .01)
  PUB <- ifelse(!is.null(P), max(P), .99)
  if(nlevels == 2){
    n0  <- ifelse(!is.null(n), mean(n), nJ0[1])
    nLB <- ifelse(!is.null(n), min(n), LB[1])
    nUB <- ifelse(!is.null(n), max(n), Inf)
    J0  <- ifelse(!is.null(J), mean(J), nJ0[2])
    JLB <- ifelse(!is.null(J), min(J), LB[2])
    JUB <- ifelse(!is.null(J), max(J), Inf)
    ssP0  <- c(n0,J0,P0)
    ssPLB <- c(nLB,JLB,PLB)
    ssPUB <- c(nUB,JUB,PUB)
  }else if(nlevels == 3){
    n0  <- ifelse(!is.null(n), mean(n), nJK0[1])
    nLB <- ifelse(!is.null(n), min(n), LB[1])
    nUB <- ifelse(!is.null(n), max(n), Inf)
    J0  <- ifelse(!is.null(J), mean(J), nJK0[2])
    JLB <- ifelse(!is.null(J), min(J), LB[2])
    JUB <- ifelse(!is.null(J), max(J), Inf)
    K0  <- ifelse(!is.null(K), mean(K), nJK0[3])
    KLB <- ifelse(!is.null(K), min(K), LB[3])
    KUB <- ifelse(!is.null(K), max(K), Inf)
    ssP0  <- c(n0,J0,K0,P0)
    ssPLB <- c(nLB,JLB,KLB,PLB)
    ssPUB <- c(nUB,JUB,KUB,PUB)
  }else if(nlevels == 4){
    n0  <- ifelse(!is.null(n), mean(n), nJKL0[1])
    nLB <- ifelse(!is.null(n), min(n), LB[1])
    nUB <- ifelse(!is.null(n), max(n), Inf)
    J0  <- ifelse(!is.null(J), mean(J), nJKL0[2])
    JLB <- ifelse(!is.null(J), min(J), LB[2])
    JUB <- ifelse(!is.null(J), max(J), Inf)
    K0  <- ifelse(!is.null(K), mean(K), nJKL0[3])
    KLB <- ifelse(!is.null(K), min(K), LB[3])
    KUB <- ifelse(!is.null(K), max(K), Inf)
    L0  <- ifelse(!is.null(L), mean(L), nJKL0[4])
    LLB <- ifelse(!is.null(L), min(L), LB[4])
    LUB <- ifelse(!is.null(L), max(L), Inf)
    ssP0  <- c(n0,J0,K0,L0,P0)
    ssPLB <- c(nLB,JLB,KLB,LLB,PLB)
    ssPUB <- c(nUB,JUB,KUB,LUB,PUB)
  }

  # constrained optimal solution
  nlopt.ssP <- nloptr::auglag(x0 = ssP0, fn = fn.min, heq = fn.constr,
                              localsolver = local.optim, localtol = 1e-8,
                              lower = ssPLB, upper = ssPUB)
  if(nlopt.ssP$convergence<0){
    warning(simpleWarning("Solution is not feasible. It cannot be trusted. Try changing default settings"))
  }
  if(all(nlopt.ssP$par == ssP0)){
    warning(simpleWarning("Solution is same as starting values. It cannot be trusted. Try changing default settings"))
  }
  if(any(nlopt.ssP$par <= 0)){
    stop(simpleError("Solution is not feasible due to non-positive sample size specifications"))
  }

  # exact solution
  ssP1 <- nlopt.ssP$par
  est.cost <- .min.cost(ssP1)
  est.mdes <- .eq.mdes(ssP1) + mdes
  est.power <- .eq.power(ssP1) + power
  exact.optim <- cbind(t(ssP1), est.cost, est.mdes, est.power)
  colnames(exact.optim) <-  c(c("n", "J", "K", "L")[1:nlevels], "P", "cost", "mdes", "power")

  # round solution
  ssP1 <- nlopt.ssP$par
  ssP1[nlevels] <- round(ssP1[nlevels])
  if(nlevels==2){
    ssP1[1] <- round(prod(ssP1[1:2]))/ssP1[2]
  }else if(nlevels==3){
    ssP1[2] <- round(prod(ssP1[2:3]))/ssP1[3]
    ssP1[1] <- round(prod(ssP1[1:3]))/prod(ssP1[2:3])
  }else if(nlevels==4){
    ssP1[3] <- round(prod(ssP1[3:4]))/ssP1[4]
    ssP1[2] <- round(prod(ssP1[2:4]))/prod(ssP1[3:4])
    ssP1[1] <- round(prod(ssP1[1:4]))/prod(ssP1[2:4])
  }
  if(nlevels==rlevel){
    ssP1[nlevels+1] <- round(ssP1[nlevels+1]*ssP1[rlevel]) / ssP1[rlevel]
  }else{
    ssP1[nlevels+1] <- round(ssP1[nlevels+1]*round(prod(ssP1[rlevel:nlevels]))) / round(prod(ssP1[rlevel:nlevels]))
  }
  est.cost <- .min.cost(ssP1)
  est.mdes <- .eq.mdes(ssP1) + mdes
  est.power <- .eq.power(ssP1) + power
  round.optim <- cbind(t(ssP1), est.cost, est.mdes, est.power)
  colnames(round.optim) <-  c(c("n", "J", "K", "L")[1:nlevels], "P", "cost", "mdes", "power")

  # grid to approximate an integer solution
  ssP1 <- nlopt.ssP$par
  lower.grid <- ssP1[1:nlevels] - rev(1:nlevels)*gm
  ssLB <- ssPLB[1:nlevels]
  ssUB <- ssPUB[1:nlevels]
  lower.grid[lower.grid <= ssLB] <- ssLB[lower.grid <= ssLB]
  upper.grid <- ssP1[1:nlevels] + rev(1:nlevels)*gm
  upper.grid[upper.grid >= ssUB] <- ssUB[upper.grid >= ssUB]
  lu.grid <- cbind(lower.grid, upper.grid)
  ss.list <- list()
  for(i in 1:nlevels){
    var <- paste("var", i, sep = "")
    ss.list[[var]] <- lu.grid[i,1]:lu.grid[i,2]
  }
  gridss <- round(as.matrix(expand.grid(ss.list)))
  gridssP <- cbind(gridss, ssP1[nlevels+1])
  if(nlevels==rlevel){
    gridssP[,nlevels+1] <- round(gridssP[,nlevels+1]*gridssP[,rlevel]) / gridssP[,rlevel]
  }else{
    gridssP[,nlevels+1] <- round(gridssP[,nlevels+1]*round(apply(gridssP[,(rlevel:nlevels)],1,prod))) / round(apply(gridssP[,(rlevel:nlevels)],1,prod))
  }
  integer.optim <- matrix(NA, nrow(gridssP), nlevels+4)
  integer.optim[,1:(nlevels+1)] <- gridssP
  for(i in 1:nrow(gridss)){
    integer.optim[i,nlevels+2] <- .min.cost(gridssP[i,])
    integer.optim[i,nlevels+3] <- .eq.mdes(gridssP[i,])
    integer.optim[i,nlevels+4] <- .eq.power(gridssP[i,])
  }
  # output result
  if(constrain == "power"){
    idx <- order(abs(integer.optim[,nlevels+4]), integer.optim[,nlevels+2],
                 decreasing = FALSE)[1:ncase]
  }else if(constrain == "mdes"){
    idx <- order(abs(integer.optim[,nlevels+3]), integer.optim[,nlevels+2],
                 decreasing = FALSE)[1:ncase]
  }else if(constrain == "cost"){
    idx <- order(abs(integer.optim[,nlevels+2]-cost), -integer.optim[,nlevels+4],
                 decreasing = FALSE)[1:ncase]
  }
  integer.optim <- integer.optim[idx,]
  integer.optim[,nlevels+3] <- integer.optim[,nlevels+3] + mdes
  integer.optim[,nlevels+4] <- integer.optim[,nlevels+4] + power
  colnames(integer.optim) <- c(c("n", "J", "K", "L")[1:nlevels], "P", "cost", "mdes", "power")

  nloptr <- list(pars = nlopt.ssP$par,
                 obj.value = nlopt.ssP$value,
                 iter = nlopt.ssP$iter,
                 global.solver = nlopt.ssP$global_solver,
                 local.solver = nlopt.ssP$local_solver,
                 convergence = nlopt.ssP$convergence,
                 message = nlopt.ssP$message)
  optim.out <- list(fun = fun,
                    parms = get("parms", parent.frame()),
                    nloptr = nloptr,
                    exact.optim = exact.optim,
                    round.optim = round.optim,
                    integer.optim = integer.optim
                    )

  return(invisible(optim.out))
}

# error handler
.error.handler <- function(x){

  # exclude default NULL arguments
  names.x <- names(x)
  idx.notnull <- match(names(lapply(x, is.null)[!lapply(x, is.null) == TRUE]), names.x)
  parms.notnull <- x[idx.notnull]
  # redefine the check list
  names.x <- names(parms.notnull)
  x <-  parms.notnull

  # validity check for sample sizes
  idx.n <- match(c("n","J","K","L"),  names.x)
  length.list.n <- length(unlist(x[idx.n[!is.na(idx.n)]]))
  length.unlist.n <- length(x[idx.n[!is.na(idx.n)]])
  if(length.list.n == length.unlist.n){
    if(any(x[idx.n[!is.na(idx.n)]] <= 0)){
      stop(
        simpleError(
          paste("Argument", sQuote(names.x[idx.n[!is.na(idx.n)]][x[idx.n[!is.na(idx.n)]] <= 0]), "cannot be negative")
        )
      )
    }
    if(any(x[idx.n[!is.na(idx.n)]] == 1)){
      warning(
        simpleWarning(
          paste("Argument", sQuote(names.x[idx.n[!is.na(idx.n)]][x[idx.n[!is.na(idx.n)]] == 1]), "is 1")
        )
      )
    }
  }

  # validity check for number of covariates
  idx.g <- match(c("g1", "g2", "g3", "g4"),  names.x)
  if(any(x[idx.g[!is.na(idx.g)]] < 0)){
    stop(
      simpleError(
        paste("Argument", sQuote(names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] < 0]), "cannot be negative")
      )
    )
  }

  # validity check for variance parameters, proportions, and probabilities
  idx.var <- match(c("R12","R22","R32","R42", "RT22", "RT32", "RT42",
                     "rho2", "rho3", "rho4", "omega2", "omega3", "omega4", "Q", "alpha", "power"),  names.x)
  if(any(x[idx.var[!is.na(idx.var)]] < 0) | any(x[idx.var[!is.na(idx.var)]] > 1)){
    stop(
      simpleError(
        paste("Argument", sQuote(names.x[idx.var[!is.na(idx.var)]][x[idx.var[!is.na(idx.var)]] < 0 | x[idx.var[!is.na(idx.var)]] > 1]), "should be between 0 and 1")
      )
    )
  }
  if(any(x$P < .01 | x$P > .99 )){
    stop(
      simpleError(
        paste("Argument", sQuote(quote(P)), "should be between .01 and .99")
      )
    )
  }

  idx.r2 <- match(c("R12","R22","R32","R42", "RT22", "RT32", "RT42"),  names.x)
  if(any(x[idx.r2[!is.na(idx.r2)]] > 0) & any(x[idx.g[!is.na(idx.g)]] == 0)){
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] > 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] == 0]
    if(substr(err.r2, 2, 2) == substr(err.g, 2, 2)){
      stop(
        simpleError(
          paste("Value of argument", sQuote(err.r2), "should correspond to the value of argument", sQuote(err.g))
        )
      )
    }
  }else if(any(x[idx.r2[!is.na(idx.r2)]] == 0) & any(x[idx.g[!is.na(idx.g)]] != 0 )){
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] > 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] == 0]
    if(substr(err.r2, 2, 2) == substr(err.g, 2, 2)){
      stop(
        simpleError(
          paste("Value of argument", sQuote(err.r2), "should correspond to the value of argument", sQuote(err.g))
        )
      )
    }
  }

  # warn for negative treatment effect
  if(any(x$mdes < 0)){
    warning(
      simpleWarning(paste("Argument", sQuote(quote(mdes)), "has negative value"))
    )
  }

  # validty check for two-tailed test
  if(x$two.tail != TRUE & x$two.tail != FALSE){
    stop(
      simpleError(paste("Value for argument", sQuote(quote(two.tail)), "should be either 'TRUE' or 1, 'FALSE'or 0"))
    )
  }

  # validty check for other cosa related parameters
  if("gm" %in% names.x){
    if(x$gm <= 0){
      stop(
        simpleError(paste("Value for argument", sQuote(quote(gm)), "should be >= 1"))
      )
    }
    if(x$ncase < 2){
      stop(
        simpleError(paste("Value for argument", sQuote(quote(ncase)), "should be >= 2"))
      )
    }
  }

}

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

compare.cosa <- function(x){
  if(inherits(x, c("parms","cosa"))){
    capture.output(invisible(suppressWarnings({
      x$parms$optimizer <- "auglag_cobyla"
      cobyla <- do.call(x$fun, x$parms)$round.optim
      x$parms$optimizer <- "auglag_mma"
      mma <- do.call(x$fun, x$parms)$round.optim
      x$parms$optimizer <- "auglag_slsqp"
      slsqp <- do.call(x$fun, x$parms)$round.optim
      x$parms$optimizer <- "auglag_lbfgs"
      lbfgs <- do.call(x$fun, x$parms)$round.optim
      cosa <- rbind(cobyla, mma, slsqp, lbfgs)
      row.names(cosa) <- c("auglag_cobyla", "auglag_mma", "auglag_slsqp", "auglag_lbfgs")
    })))
    print(round(cosa, 3))
    return(invisible(cosa))
  }else{
    stop("x should be an object returned from one of the 'cosa' functions", call.=FALSE)
  }
}
