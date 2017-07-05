mdes.cra2r2 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, P=.50, g2=0, R12=0, R22=0,
                        n, J, ...){
  df <- J - g2 - 2
  SSE <- sqrt(rho2*(1-R22)/(P*(1-P)*J) +
               (1-rho2)*(1-R12)/(P*(1-P)*J*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.cra2r2",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}
# example
# mdes.cra2r2(rho2=.20, n=4, J=20, alpha=.01)

mdes.mod1n.cra2r2 <- function(power=.80, alpha=.05, two.tail=TRUE,
                             rho2, P=.50, Q=NULL, g1=0, R12=0,
                             n, J, ...){
  df = n*J - J - g1 - 2
  SSE <- ifelse(is.null(Q),
                 sqrt((1-rho2)*(1-R12)/(P*(1-P)*J*n)), # continuous mod
                 sqrt((1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*n)) # binary mod
                 )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.mod1n.cra2r2",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(mdes.out))
}
# example
# mdes.mod1n.cra2r2(rho2=.20, n=4, J=20, Q=.3)

mdes.mod1r.cra2r2 <- function(power=.80, alpha=.05, two.tail=TRUE,
                              rho2, omega2, RT22=0, P=.50, Q=NULL, g1=0, R12=0,
                              n, J, ...){
  df = J - g1 - 2
  SSE <- ifelse(is.null(Q),
                 sqrt(rho2*omega2*(1-RT22)/(P*(1-P)*J) +
                    (1-rho2)*(1-R12)/(P*(1-P)*J*n)), # continuous mod
                 sqrt(rho2*omega2*(1-RT22)/(P*(1-P)*J) +
                    (1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*n)) # binary mod
                 )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.mod1r.cra2r2",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(mdes.out))
}
# example
# mdes.mod1r.cra2r2(rho2=.2, omega2=.2, RT22=.2, n=4, J=20)


mdes.mod2.cra2r2 <- function(power=.80, alpha=.05, two.tail=TRUE,
                             rho2, P=.50, Q=NULL, g2=0, R12=0, R22=0,
                             n, J, ...){
  df = J - g2 - 4
  SSE <- ifelse(is.null(Q),
                sqrt(rho2*(1-R22)/(P*(1-P)*(J-g2-4)) +
                       (1-rho2)*(1-R12)/(P*(1-P)*(J-g2-4)*n)), # continuous mod
                sqrt(rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*(J-g2-4)) +
                       (1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*(J-g2-4)*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.mod2.cra2r2",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(mdes.out))
}
# example
# mdes.mod2.cra2r2(rho2=.20, n=4, J=20, Q=.5)


power.cra2r2 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, g2=0, P=.50, R12=0, R22=0,
                         n, J, ...){
  df <- J - g2 - 2
  SSE <- sqrt(rho2*(1-R22)/(P*(1-P)*J) +
                    (1-rho2)*(1-R12)/(P*(1-P)*J*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.cra2r2",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.cra2r2(rho2=.20, n=4, J=20)

power.mod1n.cra2r2 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                              rho2, P=.50, Q=NULL, g1=0, R12=0,
                              n, J, ...){
  df = n*J - J - g1 - 2
  SSE <- ifelse(is.null(Q),
                sqrt((1-rho2)*(1-R12)/(P*(1-P)*J*n)), # continuous mod
                sqrt((1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.mod1n.cra2r2",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(power.out))
}
# example
# power.mod1n.cra2r2(rho2=.20, n=4, J=20)

power.mod1r.cra2r2 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                              rho2, omega2, RT22=0, P=.50, Q=NULL, g1=0, R12=0,
                              n, J, ...){
  df = J - g1 - 2
  SSE <- ifelse(is.null(Q),
                sqrt(rho2*omega2*(1-RT22)/(P*(1-P)*J) +
                       (1-rho2)*(1-R12)/(P*(1-P)*J*n)), # continuous mod
                sqrt(rho2*omega2*(1-RT22)/(P*(1-P)*J) +
                       (1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.mod1r.cra2r2",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(power.out))
}
# example
# power.mod1r.cra2r2(rho2=.2, omega2=.2, RT22=.2, n=4, J=20)


power.mod2.cra2r2 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                             rho2, P=.50, Q=NULL, g2=0, R12=0, R22=0,
                             n, J, ...){
  df = J - g2 - 4
  SSE <- ifelse(is.null(Q),
                sqrt(rho2*(1-R22)/(P*(1-P)*(J-g2-4)) +
                       (1-rho2)*(1-R12)/(P*(1-P)*(J-g2-4)*n)), # continuous mod
                sqrt(rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*(J-g2-4)) +
                       (1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*(J-g2-4)*n)) # binary mod
  )
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.mod2.cra2r2",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  print(ifelse(is.null(Q),"Note: Continuous moderator","Note: Binary moderator"))
  return(invisible(power.out))
}
# example
# power.mod2.cra2r2(rho2=.20, n=4, J=20)


optimal.cra2r2 <- function(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                                power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                nJ0=c(10,10), ncase=10, gm=2,
                                constrain="power", optimizer="auglag_cobyla",
                                rho2, P=.50, g2=0, R12=0, R22=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.cra2r2"
  LB <- c(1, g2+2+1)
  df <- quote(J - g2 - 2)
  SSE <- quote(sqrt(rho2*(1-R22)/(P*(1-P)*J) +
                      (1-rho2)*(1-R12)/(P*(1-P)*J*n)))
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# optimal.cra2r2(cn=1, cJ=10, cost=560, constrain="cost", rho2=.20)

optimal.mod1n.cra2r2 <- function(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                 nJ0=c(10,10), ncase=10, gm=2,
                                 constrain="power", optimizer="auglag_cobyla",
                                 rho2, P=.50, Q=NULL, g1=0, R12=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.mod1n.cra2r2"
  LB <- c(2, g1+2+1)
  df <- quote(n*J - J - g1 - 2)
  SSE <- if(is.null(Q)){
    quote(sqrt((1-rho2)*(1-R12)/(P*(1-P)*J*n))) # continuous mod
  }else{
    quote(sqrt((1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*n))) # binary mod
  }
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# optimal.mod1n.cra2r2(cn=1, cJ=10, cost=560, constrain="cost", rho2=.20)

optimal.mod1r.cra2r2 <- function(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                                 power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                 nJ0=c(10,10), ncase=10, gm=2,
                                 constrain="power", optimizer="auglag_cobyla",
                                 rho2, omega2, P=.50, Q=NULL, g1=0, R12=0, RT22=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.mod1r.cra2r2"
  LB <- c(1, g1+4+1)
  df <- quote(J - g1 - 4)
  SSE <- if(is.null(Q)){
    quote(sqrt(rho2*omega2*(1-RT22)/(P*(1-P)*J) +
                 (1-rho2)*(1-R12)/(P*(1-P)*J*n))) # continuous mod
  }else{
    quote(sqrt(rho2*omega2*(1-RT22)/(P*(1-P)*J) +
                 (1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*J*n))) # binary mod
  }
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# optimal.mod1r.cra2r2(cn=1, cJ=10, cost=560, constrain="cost", rho2=.20, omega2=.20)


optimal.mod2.cra2r2 <- function(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                                power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                                nJ0=c(10,10), ncase=10, gm=2,
                                constrain="power", optimizer="auglag_cobyla",
                                rho2, P=.50, Q=NULL, g2=0, R12=0, R22=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.mod2.cra2r2"
  LB <- c(1, g2+4+1)
  df <- quote(J - g2 - 4)
  SSE <- if(is.null(Q)){
    quote(sqrt(rho2*(1-R22)/(P*(1-P)*(J-g2-4)) +
                 (1-rho2)*(1-R12)/(P*(1-P)*(J-g2-4)*n))) # continuous mod
  }else{
    quote(sqrt(rho2*(1-R22)/(P*(1-P)*Q*(1-Q)*(J-g2-4)) +
                 (1-rho2)*(1-R12)/(P*(1-P)*Q*(1-Q)*(J-g2-4)*n))) # binary mod
  }
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}
# example
# optimal.mod2.cra2r2(cn=1, cJ=10, cost=560, constrain="cost", rho2=.20)