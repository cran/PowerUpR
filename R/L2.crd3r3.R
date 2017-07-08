mdes.crd3r3 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, rho3, P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                        g3=0, R12=0, R22=0, R32=0, n, J, K, ...){
  df <- K - g3 - 2
  D <- .D.fun(P, k1, k2, dist.Z, RTZ)
  SSE<- sqrt(D*(rho3*(1-R32)/(P*(1-P)*K) +
                rho2*(1-R22)/(P*(1-P)*J*K) +
               (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n)))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.crd3r3",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.crd3r3(rho3=.06, rho2=.17, n=15, J=3, K=60)


power.crd3r3 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, rho3, P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                         g3=0, R12=0, R22=0, R32=0, n, J, K, ...){
  df <- K - g3 - 2
  D <- .D.fun(P, k1, k2, dist.Z, RTZ)
  SSE <- sqrt(D*(rho3*(1-R32)/(P*(1-P)*K) +
                     rho2*(1-R22)/(P*(1-P)*J*K) +
                    (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n)))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.crd3r3",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.crd3r3(mdes=0.447, rho3=.06, rho2=.17, n=15, J=3, K=60)

cosa.crd3r3 <- function(cn=0, cJ=0, cK=0, cost=NULL,
                           n=NULL, J=NULL, K=NULL, nJK0=c(10,10,10),
                           constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           rho2, rho3,
                           P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                           g3=0, R12=0, R22=0, R32=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.crd3r3"
  LB <- c(1, 1, g3+2+1)
  df <- quote(K - g3 - 2)
  SSE <- quote(sqrt(D*(rho3*(1-R32)/(P*(1-P)*K) +
                     rho2*(1-R22)/(P*(1-P)*J*K) +
                    (1-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*n))))
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

# examples

# cosa given total cost
# cosa.crd3r3(mdes=0.447, cn=1, cJ=10, cK=100, cost=5600, constrain="cost", rho2=.20, rho3=.10)

# cosa given per unit costs and power
# cosa.crd3r3(mdes=0.447, cn=1, cJ=10, cK=100, power=.80, constrain="power", rho2=.20, rho3=.10)

# cosa given per unit costs and mdes
# cosa.crd3r3(mdes=0.447, cn=1, cJ=10, cK=100, power=.80, constrain="mdes",  rho2=.20, rho3=.10)

# cosa (fixed n=10) given total cost
# cosa.crd3r3(mdes=0.447, cn=1, cJ=10, cK=100, n=10, cost=5600, constrain="cost",  rho2=.20, rho3=.10)

# cosa (fixed n=10) given per unit costs and power
# cosa.crd3r3(mdes=0.447, cn=1, cJ=10, cK=100, n=10, constrain="power", rho2=.20, rho3=.10)

# cosa (fixed n=10) given per unit costs and mdes
# cosa.crd3r3(mdes=0.447, cn=1, cJ=10, cK=100, constrain="mdes", rho3=.06, rho2=.17, n=15, J=3, K=60)

