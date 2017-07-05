mdes.bcra3r2 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, rho3, omega3, P=.50, g3=0, R12=0, R22=0, RT32=0,
                        n, J, K, ...){
  df <- K - g3 - 1
  SSE <- sqrt(rho3*omega3*(1-RT32)/K +
                rho2*(1-R22)/(P*(1-P)*J*K) +
               (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.bcra3r2",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.bcra3r2(rho3=.13, rho2=.10, omega3=.4, n=10, J=6, K=24)


power.bcra3r2 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, rho3,  omega3, P=.50, g3=0, R12=0, R22=0, RT32=0,
                         n, J, K, ...){
  df <- K - g3 - 1
  SSE <- sqrt(rho3*omega3*(1-RT32)/K +
                     rho2*(1-R22)/(P*(1-P)*J*K) +
                    (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.bcra3r2",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.bcra3r2(rho3=.13, rho2=.10, omega3=.4, n=10, J=6, K=24)

cosa.bcra3r2 <- function(cn=0, cJ=0, cK=0, cost=NULL,
                           n=NULL, J=NULL, K=NULL, P=NULL,
                           nJK0=c(10,10,10), P0=.50,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                           rho2, rho3, omega3, g3=0, R12=0, R22=0, RT32=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.bcra3r2"
  LB <- c(1, 1, g3+1+1)
  df <- quote(K - g3 - 1)
  SSE <- quote(sqrt(rho3*omega3*(1-RT32)/K +
                rho2*(1-R22)/(P*(1-P)*J*K) +
                (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n)))
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

optimal.bcra3r2 <- function(...){
  .Deprecated(new="cosa.bcra3r2")
  cosa.bcra3r2(...)
}

# examples

# cosa given total cost
# cosa.bcra3r2(cn=1, cJ=10, cK=100, cost=5600, constrain="cost", rho3=.15, rho2=.20, omega3=.50)

# cosa given per unit costs and power
# cosa.bcra3r2(cn=1, cJ=10, cK=100, power=.80, constrain="power", rho3=.15, rho2=.20, omega3=.50)

# cosa given per unit costs and mdes
# cosa.bcra3r2(cn=1, cJ=10, cK=100, power=.80, constrain="mdes", rho3=.15, rho2=.20, omega3=.50)

# cosa (fixed n=10) given total cost
# cosa.bcra3r2(cn=1, cJ=10, cK=100, n=10, cost=5600, constrain="cost", rho3=.15, rho2=.20, omega3=.50)

# cosa (fixed n=10) given per unit costs and power
# cosa.bcra3r2(cn=1, cJ=10, cK=100, n=10, constrain="power", rho3=.15, rho2=.20, omega3=.50)

# cosa (fixed n=10) given per unit costs and mdes
# cosa.bcra3r2(cn=1, cJ=10, cK=100, n=10, constrain="mdes", rho3=.15, rho2=.20, omega3=.50)


