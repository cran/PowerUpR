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
                      ncase = 10, gm = 2, constrain = "power", optimizer = "auglag_slsq",
                      rho2, rho3, rho4, omega2, omega3, omega4,
                      Q=NULL, g1=0, g2=0, g3=0, g4=0,
                      R12=0, R22=0, R32=0, R42=0, RT22=0, RT32=0, RT42=0, R2T2=0, R3T2=0,
                      ...){

  df <- get("df", parent.frame())
  SSE <- get("SSE", parent.frame())
  two.tail <- get("two.tail", parent.frame())
  alpha <- get("alpha", parent.frame())
  power <- get("power", parent.frame())
  mdes <- get("mdes", parent.frame())
  LB <- get("LB", parent.frame())
  fun <- get("fun", parent.frame())

  fun.parsed <- scan(text = fun, what = "character", sep=".", quiet = TRUE)
  rlevel <- as.numeric(substr(fun, nchar(fun), nchar(fun)))
  nlevels <- as.numeric(substr(fun, nchar(fun)-2, nchar(fun)-2)) # nlevels <- length(LB)

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

  cost.list <- list(cn, cJ, cK, cL)
  cost.names <- c("cn", "cJ", "cK", "cL")
  cost.lengths <- unlist(lapply(cost.list, length))
  if(rlevel < nlevels){
    if(any(cost.lengths[rlevel:nlevels] == 2)){
      stop("Unequal cost applies to levels at or below randomization (or discontinuity) level", call.=FALSE)
    }else if(any(cost.lengths[1:rlevel] > 2)){
      stop("Marginal costs cannot have a length greater than two", call.=FALSE)
    }
  }else{
    if(any(cost.lengths[1:nlevels] == 2)){
      message("The first value in marginal costs refers to per unit in treatment condition")
    }else if(any(cost.lengths[1:nlevels] > 2)){
      stop("Marginal costs cannot have a length greater than two", call.=FALSE)
    }
  }

  if(!is.null(cost)){
    if(length(cost)>1){
      stop("Total cost cannot have a length greater than one", call.=FALSE)
    }else if(cost<0){
      stop(" Total cost cannot have a negative value", call.=FALSE)
    }
  }

  if(fun.parsed[2] %in% c("crd2r2", "crd3r3", "bird2r1", "bird2f1", "bcrd3r2", "bcrd3f2")){
    if(length(P)>1 | !is.null(P)){
      stop("P cannot be NULL or have a length greater than one in RDD", call.=FALSE)
    }
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

  if(length(cn)==1){
    cn <- c(cn,cn)
  }
  if(length(cJ)==1){
    cJ <- c(cJ,cJ)
  }
  if(length(cK)==1){
    cK <- c(cK,cK)
  }
  if(length(cL)==1){
    cL <- c(cL,cL)
  }

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

  # constrained optimal sample allocation
  nlopt.ssP <- nloptr::auglag(x0 = ssP0, fn = fn.min, heq = fn.constr,
                              localsolver = local.optim, localtol = 1e-8,
                              lower = ssPLB, upper = ssPUB)
  if(nlopt.ssP$convergence < 0 | all(nlopt.ssP$par == ssP0) | any(nlopt.ssP$par <= 0)){
    stop("Solution is not feasible. Change default settings", call.=FALSE)
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

  # integer approximations
  ssP1 <- nlopt.ssP$par
  ssLB <- ssPLB[1:nlevels]
  ssUB <- ssPUB[1:nlevels]
  lower.grid <- ssP1[1:nlevels] - rev(1:nlevels)*gm
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

