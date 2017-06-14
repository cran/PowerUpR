mdes.bcra4f3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, rho3, P=.50, R12=0, R22=0, R32=0, g3=0,
                        n, J, K, L, ...){
  df <- L * (K - 2) - g3
  SSE <- sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                  rho2*(1-R22)/(P*(1-P)*J*K*L) +
                  (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.bcra4f3",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.bcra4f3(alpha=.05, two.tail=TRUE, power=.80, rho3=.15, rho2=.15, n=10, J=4, K=4, L=15)

power.bcra4f3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, rho3, P=.50, R12=0, R22=0, R32=0, g3=0,
                         n, J, K, L, ...){
  df <- L * (K - 2) - g3
  SSE <- sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                     rho2*(1-R22)/(P*(1-P)*J*K*L) +
                     (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.bcra4f3",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.bcra4f3(mdes=0.24, alpha=.05, two.tail=TRUE, rho3=.15, rho2=.15, n=10, J=4, K=4, L=15)

optimal.bcra4f3 <- function(cn=0, cJ=0, cK=0, cL=0, cost=NULL, n=NULL, J=NULL, K=NULL, L=NULL,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           nJKL0=c(10,10,10,10), ncase=10, gm=2,
                           constrain="power", optimizer="auglag_cobyla",
                           rho3, rho2, P=.50, g3=0,  R32=0, R22=0, R12=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.bcra4f3"
  LB <- c(1, 1, 3, g3+1)
  df <- quote(L * (K - 2) - g3)
  SSE <- quote(sqrt(rho3*(1-R32)/(P*(1-P)*K*L) +
                rho2*(1-R22)/(P*(1-P)*J*K*L) +
                (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n)))
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
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
