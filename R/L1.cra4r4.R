mdes.cra4r4 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, rho3, rho4, P=.50, R12=0, R22=0, R32=0, R42=0, g4=0,
                        n, J, K, L, ...){
  df <- L - g4 - 2
  SSE <- sqrt(rho4*(1-R42)/(P*(1-P)*L) +
                rho3*(1-R32)/(P*(1-P)*K*L) +
                rho2*(1-R22)/(P*(1-P)*J*K*L) +
               (1-rho4-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.cra4r4",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.cra4r4(rho4=.05, rho3=.05, rho2=.10, n=10, J=2, K=3, L=20)

power.cra4r4 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, rho3, rho4, P=.50, R12=0, R22=0, R32=0, R42=0, g4=0,
                         n, J, K, L, ...){
  df <- L - g4 - 2
  SSE <- sqrt(rho4*(1-R42)/(P*(1-P)*L) +
                     rho3*(1-R32)/(P*(1-P)*K*L) +
                     rho2*(1-R22)/(P*(1-P)*J*K*L) +
                    (1-rho4-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.cra4r4",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.cra4r4(rho4=.05, rho3=.05, rho2=.10, n=10, J=2, K=3, L=20)


cosa.cra4r4 <- function(cn=0, cJ=0, cK=0, cL=0, cost=NULL,
                           n=NULL, J=NULL, K=NULL, L=NULL, P=NULL,
                           nJKL0=c(10,10,10,10), P0=.50,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                           rho4, rho3, rho2, g4=0, R42=0, R32=0, R22=0, R12=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.cra4r4"
  LB <- c(1, 1, 1, g4+2+1)
  df <- quote(L - g4 - 2)
  SSE <- quote(sqrt(rho4*(1-R42)/(P*(1-P)*L) +
                rho3*(1-R32)/(P*(1-P)*K*L) +
                rho2*(1-R22)/(P*(1-P)*J*K*L) +
                (1-rho4-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n)))
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

optimal.cra4r4 <- function(...){
  .Deprecated(new="cosa.cra4r4")
  cosa.cra4r4(...)
}

# examples

# cosa given total cost
# cosa.cra4r4(cn=1, cJ=10, cK=100, cL=1000, cost=75600, constrain="cost", rho4=.10, rho3=.15, rho2=.20)

# cosa given per unit costs and power
# cosa.cra4r4(cn=1, cJ=10, cK=100, cL=1000, power=.80, constrain="power", rho4=.10, rho3=.15, rho2=.20)

# cosa given per unit costs and mdes
# cosa.cra4r4(cn=1, cJ=10, cK=100, cL=1000, power=.80, constrain="mdes", rho4=.10, rho3=.15, rho2=.20)

# cosa (fixed n=10, and J=2) given total cost
# cosa.cra4r4(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, cost=75600, constrain="cost", rho4=.10, rho3=.15, rho2=.20)

# cosa (fixed n=10, and J=2) given per unit costs and power
# cosa.cra4r4(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, constrain="power", rho4=.10, rho3=.15, rho2=.20, nJKL0=c(10,2,3,3))

# cosa (fixed n=10, and J=2) given per unit costs and mdes
# cosa.cra4r4(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, constrain="mdes", rho4=.10, rho3=.15, rho2=.20)
