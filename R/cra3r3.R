mdes.cra3r3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                             rho2, rho3, P=.50, g3=0, R12=0, R22=0, R32=0,
                             n, J, K, ...){
  df <- K - g3 - 2
  SSE <- sqrt(rho3*(1-R32)/(P*(1-P)*K) +
                rho2*(1-R22)/(P*(1-P)*J*K) +
                (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.cra3r3",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.cra3r3(rho3=.06, rho2=.17, n=15, J=3, K=60)

mdes.mod1n.cra3r3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                             rho2, rho3, P=.50, Q=NULL, g1=0, R12=0,
                             n, J, K, ...){
  df <- n*J*K - J*K - K - g1 - 2
  SSE <- ifelse(is.null(Q),
                sqrt((1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n)), # continuous mod
                sqrt((1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*K*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.mod1n.cra3r3",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(mdes.out))
}

# example
# mdes.mod1n.cra3r3(rho3=.1, rho2=.1, Q=.5, g1=1, R12=.3, n=20, J=4, K=60)
# mdes.mod1n.cra3r3(rho3=.1, rho2=.1, g1=1, R12=.3, n=20, J=4, K=60)


mdes.mod1r.cra3r3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                              rho2, rho3, omega2, omega3, P=.50, Q=NULL, g1=0, R12=0, RT22=0, RT32=0,
                              n, J, K, ...){
  df <- K - g1 - 1
  SSE <- ifelse(is.null(Q),
                sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                       rho2*omega2*(1-RT22)/(P*(1-P)*K*J) +
                       (1-rho2-rho3)*(1-R12)/(P*(1-P)*K*J*n)), # continuous mod
                sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                       rho2*omega2*(1-RT22)/(P*(1-P)*K*J) +
                       (1-rho2-rho3)*(1-R12)/(P*(1-P)*Q*(1-Q)*K*J*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.mod1r.cra3r3",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(mdes.out))
}

# example
# mdes.mod1r.cra3r3(rho3=.05, rho2=.12, omega2=.08, omega3=.07, P=.4, Q=.7, g1=1, R12=.20, RT22=0, RT32=0,  n=20, J=4, K=60)
# mdes.mod1r.cra3r3(rho3=.05, rho2=.12, omega2=.08, omega3=.07, P=.4, g1=1, R12=.20, RT22=0, RT32=0,  n=20, J=4, K=60)


mdes.mod2n.cra3r3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                              rho2, rho3, P=.50, Q=NULL, g2=0, R12=0, R22=0,
                              n, J, K, ...){
  df <- J*K - K - g2 - 2
  SSE <- ifelse(is.null(Q),
                sqrt(rho2*(1-R22)/(P*(1-P)*J*K) +
                      (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n)), # continuous mod
                sqrt(rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*J*K) +
                      (1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*K*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.mod2n.cra3r3",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(mdes.out))
}

# example
# mdes.mod2n.cra3r3(rho3=.1, rho2=.1, Q=.5, g2=1, R12=.3, R22=.4,  n=20, J=4, K=60)
# mdes.mod2n.cra3r3(rho3=.1, rho2=.1, g2=1, R12=.3, R22=.4, n=20, J=4, K=60)


mdes.mod2r.cra3r3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                              rho2, rho3, omega3, P=.50, Q=NULL, g2=0, R12=0, R22=0, RT32=0,
                              n, J, K, ...){
  df <- K - g2 - 1
  SSE <- ifelse(is.null(Q),
                sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                       rho2*(1-R22)/(P*(1-P)*K*J) +
                       (1-rho2-rho3)*(1-R12)/(P*(1-P)*K*J*n)), # continuous mod
                sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                       rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*K*J) +
                       (1-rho2-rho3)*(1-R12)/(P*(1-P)*Q*(1-Q)*K*J*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.mod2r.cra3r3",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(mdes.out))
}

# example
# mdes.mod2r.cra3r3(rho3=.1, rho2=.1, omega3=.05, Q=.5, g2=1, R12=.30, R22=.4, RT32=0,  n=20, J=4, K=60)
# mdes.mod2r.cra3r3(rho3=.1, rho2=.1, omega3=.05, g2=1, R12=.30, R22=.4, RT32=0,  n=20, J=4, K=60)

mdes.mod3.cra3r3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                             rho2, rho3, P=.50, Q=NULL, g3=0, R12=0, R22=0, R32=0,
                             n, J, K, ...){
  df <- K-g3-4
  SSE <- ifelse(is.null(Q),
                sqrt(rho3*(1-R32)/(P*(1-P)*(K-g3-4)) +
                       rho2*(1-R22)/(P*(1-P)*J*(K-g3-4)) +
                       (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*(K-g3-4)*n)), # continuous mod
                sqrt(rho3*(1-R32)/(P*(1-P)*Q*(1-Q)*(K-g3-4)) +
                       rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*J*(K-g3-4)) +
                       (1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*(K-g3-4)*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.mod3.cra3r3",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(mdes.out))
}

# example
# mdes.mod3.cra3r3(rho3=.1, rho2=.1, omega3=.05, Q=.5, g3=1, R12=.3, R22=.4, R32=.5, n=20, J=4, K=60)
# mdes.mod3.cra3r3(rho3=.1, rho2=.1, omega3=.05, g3=1, R12=.3, R22=.4, R32=.5, n=20, J=4, K=60)

power.cra3r3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                             rho2, rho3, P=.50, g3=0, R12=0, R22=0, R32=0,
                             n, J, K, ...){
  df <- K - g3 - 2
  SSE <- sqrt(rho3*(1-R32)/(P*(1-P)*K) +
                rho2*(1-R22)/(P*(1-P)*J*K) +
                (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.cra3r3",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}

# example
# power.cra3r3(mdes=.269, rho3=.06, rho2=.17, n=15, J=3, K=60)

power.mod1n.cra3r3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                              rho2, rho3, P=.50, Q=NULL, g1=0, R12=0,
                              n, J, K, ...){
  df <- n*J*K - J*K - K - g1 - 2
  SSE <- ifelse(is.null(Q),
                sqrt((1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n)), # continuous mod
                sqrt((1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*K*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  .error.handler(parms)
  power.out <-  list(fun = "power.mod1n.cra3r3",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(power.out))
}

# example
# power.mod1n.cra3r3(mdes=.12, rho3=.1, rho2=.1, Q=.5, g1=1, R12=.3, n=20, J=4, K=60)
# power.mod1n.cra3r3(mdes=.06, rho3=.1, rho2=.1, g1=1, R12=.3, n=20, J=4, K=60)


power.mod1r.cra3r3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                              rho2, rho3, omega2, omega3, P=.50, Q=NULL, g1=0, R12=0, RT22=0, RT32=0,
                              n, J, K, ...){
  df <- K - g1 - 1
  SSE <- ifelse(is.null(Q),
                sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                       rho2*omega2*(1-RT22)/(P*(1-P)*K*J) +
                       (1-rho2-rho3)*(1-R12)/(P*(1-P)*K*J*n)), # continuous mod
                sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                       rho2*omega2*(1-RT22)/(P*(1-P)*K*J) +
                       (1-rho2-rho3)*(1-R12)/(P*(1-P)*Q*(1-Q)*K*J*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.mod1r.cra3r3",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(power.out))
}

# example
# power.mod1r.cra3r3(mdes=.16, rho3=.05, rho2=.12, omega2=.08, omega3=.07, P=.4, Q=.7, g1=1, R12=.20, RT22=0, RT32=0,  n=20, J=4, K=60)
# power.mod1r.cra3r3(mdes=.09, rho3=.05, rho2=.12, omega2=.08, omega3=.07, P=.4, g1=1, R12=.20, RT22=0, RT32=0,  n=20, J=4, K=60)


power.mod2n.cra3r3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                              rho2, rho3, P=.50, Q=NULL, g2=0, R12=0, R22=0,
                              n, J, K, ...){
  df <- J*K - K - g2 - 2
  SSE <- ifelse(is.null(Q),
                sqrt(rho2*(1-R22)/(P*(1-P)*J*K) +
                       (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n)), # continuous mod
                sqrt(rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*J*K) +
                       (1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*K*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  .error.handler(parms)
  power.out <-  list(fun = "power.mod2n.cra3r3",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(power.out))
}

# example
# power.mod2n.cra3r3(mdes=.22, rho3=.1, rho2=.1, Q=.5, g1=1, R12=.3, R22=.4,  n=20, J=4, K=60)
# power.mod2n.cra3r3(mdes=.11, rho3=.1, rho2=.1, g1=1, R12=.3, R22=.4, n=20, J=4, K=60)


power.mod2r.cra3r3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                              rho2, rho3, omega3, P=.50, Q=NULL, g2=0, R12=0, R22=0, RT32=0,
                              n, J, K, ...){
  df <- K - g2 - 1
  SSE <- ifelse(is.null(Q),
                sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                       rho2*(1-R22)/(P*(1-P)*K*J) +
                       (1-rho2-rho3)*(1-R12)/(P*(1-P)*K*J*n)), # continuous mod
                sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                       rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*K*J) +
                       (1-rho2-rho3)*(1-R12)/(P*(1-P)*Q*(1-Q)*K*J*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.mod2r.cra3r3",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(power.out))
}

# example
# power.mod2r.cra3r3(mdes=.22, rho3=.1, rho2=.1, omega3=.05, Q=.5, g1=1, R12=.30, R22=.4, RT32=0,  n=20, J=4, K=60)
# power.mod2r.cra3r3(mdes=.11, rho3=.1, rho2=.1, omega3=.05, g1=1, R12=.30, R22=.4, RT32=0,  n=20, J=4, K=60)

power.mod3.cra3r3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                             rho2, rho3, P=.50, Q=NULL, g3=0, R12=0, R22=0, R32=0,
                             n, J, K, ...){
  df <- K-g3-4
  SSE <- ifelse(is.null(Q),
                sqrt(rho3*(1-R32)/(P*(1-P)*(K-g3-4)) +
                       rho2*(1-R22)/(P*(1-P)*J*(K-g3-4)) +
                       (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*(K-g3-4)*n)), # continuous mod
                sqrt(rho3*(1-R32)/(P*(1-P)*Q*(1-Q)*(K-g3-4)) +
                       rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*J*(K-g3-4)) +
                       (1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*(K-g3-4)*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.mod3.cra3r3",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(power.out))
}

# example
# power.mod3.cra3r3(mdes=.30, rho3=.1, rho2=.1, omega3=.05, Q=.5, g3=1, R12=.3, R22=.4, R32=.5, n=20, J=4, K=60)
# power.mod3.cra3r3(mdes=.15, rho3=.1, rho2=.1, omega3=.05, g3=1, R12=.3, R22=.4, R32=.5, n=20, J=4, K=60)

cosa.cra3r3 <- function(cn=0, cJ=0, cK=0, cost=NULL,
                                n=NULL, J=NULL, K=NULL, P=NULL,
                                nJK0=c(10,10,10), P0=.50,
                                power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                                rho2, rho3, g3=0, R12=0, R22=0, R32=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.cra3r3"
  LB <- c(1, 1, g3+2+1)
  df <- quote(K - g3 - 2)
  SSE <- quote(sqrt(rho3*(1-R32)/(P*(1-P)*K) +
                rho2*(1-R22)/(P*(1-P)*J*K) +
                (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n)))
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

optimal.cra3r3 <- function(...){
  .Deprecated(new="cosa.cra3r3")
  cosa.cra3r3(...)
}

# example
# cosa.cra3r3(cn=1, cJ=10, cK=100, n=10, constrain="mdes", rho2=.20, rho3=.10)



cosa.mod1n.cra3r3 <- function(cn=0, cJ=0, cK=0, cost=NULL,
                                 n=NULL, J=NULL, K=NULL, P=NULL,
                                 nJK0=c(10,10,10), P0=.50,
                                 constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                 rho2, rho3, Q=NULL, g1=0, R12=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.mod1n.cra3r3"
  LB <- c(2, 2, g1+2+1)
  df <- quote(n*J*K - J*K - K - g1 - 2)
  SSE <- if(is.null(Q)){
    quote(sqrt((1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n))) # continuous mod
  }else{
    quote(sqrt((1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*K*n))) # binary mod
  }
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# cosa.mod1n.cra3r3(cn=1, cJ=10, cK=100, n=10, constrain="mdes", rho2=.20, rho3=.10)

cosa.mod1r.cra3r3 <- function(cn=0, cJ=0, cK=0, cost=NULL,
                                 n=NULL, J=NULL, K=NULL, P=NULL,
                                 nJK0=c(10,10,10), P0=.50,
                                 constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                 rho2, rho3, omega2, omega3, Q=NULL, g1=0, R12=0, RT22=0, RT32=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.mod1r.cra3r3"
  LB <- c(1, 1, g1+1+1)
  df <- quote(K - g1 - 1)
  SSE <- if(is.null(Q)){
    quote(sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                 rho2*omega2*(1-RT22)/(P*(1-P)*K*J) +
                 (1-rho2-rho3)*(1-R12)/(P*(1-P)*K*J*n))) # continuous mod
  }else{
    quote(sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                 rho2*omega2*(1-RT22)/(P*(1-P)*K*J) +
                 (1-rho2-rho3)*(1-R12)/(P*(1-P)*Q*(1-Q)*K*J*n))) # binary mod
  }
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# cosa.mod1r.cra3r3(cn=1, cJ=10, cK=100, n=10, constrain="mdes", rho2=.20, rho3=.10, omega2=.10, omega3=.10)

cosa.mod2n.cra3r3 <- function(cn=0, cJ=0, cK=0, cost=NULL,
                                 n=NULL, J=NULL, K=NULL, P=NULL,
                                 nJK0=c(10,10,10), P0=.50,
                                 constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                 rho2, rho3, Q=NULL, g2=0, R12=0, R22=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.mod2n.cra3r3"
  LB <- c(1, 2, g2+2+1)
  df <- quote(J*K - K - g2 - 2)
  SSE <- if(is.null(Q)){
    quote(sqrt(rho2*(1-R22)/(P*(1-P)*J*K) +
                 (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n))) # continuous mod
  }else{
    quote(sqrt(rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*J*K) +
                 (1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*K*n))) # binary mod
  }
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# cosa.mod2n.cra3r3(cn=1, cJ=10, cK=100, n=10, constrain="mdes", rho2=.20, rho3=.10)

cosa.mod2r.cra3r3 <- function(cn=0, cJ=0, cK=0, cost=NULL,
                                 n=NULL, J=NULL, K=NULL, P=NULL,
                                 nJK0=c(10,10,10), P0=.50,
                                 constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                 rho2, rho3, omega3, Q=NULL, g2=0, R12=0, R22=0, RT32=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.mod2r.cra3r3"
  LB <- c(1, 1, g2+1+1)
  df <- quote(K - g2 - 1)
  SSE <- if(is.null(Q)){
    quote(sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                 rho2*(1-R22)/(P*(1-P)*K*J) +
                 (1-rho2-rho3)*(1-R12)/(P*(1-P)*K*J*n))) # continuous mod
  }else{
    quote(sqrt(rho3*omega3*(1-RT32)/(P*(1-P)*K) +
                 rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*K*J) +
                 (1-rho2-rho3)*(1-R12)/(P*(1-P)*Q*(1-Q)*K*J*n))) # binary mod
  }
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# cosa.mod2r.cra3r3(cn=1, cJ=10, cK=100, constrain="mdes", rho2=.20, rho3=.10, omega3=.10)

cosa.mod3.cra3r3 <- function(cn=0, cJ=0, cK=0, cost=NULL,
                                n=NULL, J=NULL, K=NULL, P=NULL,
                                nJK0=c(10,10,10), P0=.50,
                                constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                                power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                rho2, rho3, Q=NULL, g3=0, R12=0, R22=0, R32=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.mod3.cra3r3"
  LB <- c(1, 1, g3+4+1)
  df <- quote(K - g3 - 4)
  SSE <- if(is.null(Q)){
    quote(sqrt(rho3*(1-R32)/(P*(1-P)*(K - g3 - 4)) +
                 rho2*(1-R22)/(P*(1-P)*J*(K - g3 - 4)) +
                 (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*(K - g3 - 4)*n))) # continuous mod
  }else{
    quote(sqrt(rho3*(1-R32)/(P*(1-P)*Q*(1-Q)*(K - g3 - 4)) +
                 rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*J*(K - g3 - 4)) +
                 (1-rho3-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*(K - g3 - 4)*n))) # binary mod
  }
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# cosa.mod3.cra3r3(cn=1, cJ=10, cK=100, n=10, constrain="mdes", rho2=.20, rho3=.10)
