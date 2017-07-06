mdes.crd2r2 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                        g2=0, R12=0, R22=0, n, J, ...){
  df <- J - g2 - 2
  D <- .D.fun(P, k1, k2, dist.Z, RTZ)
  SSE <- sqrt(D*(rho2*(1-R22)/(P*(1-P)*J) +
               (1-rho2)*(1-R12)/(P*(1-P)*J*n)))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.crd2r2",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.crd2r2(rho2=.20, n=4, J=20)


power.crd2r2 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, g2=0, P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                         R12=0, R22=0, n, J, ...){
  df <- J - g2 - 2
  D <- .D.fun(P, k1, k2, dist.Z, RTZ)
  SSE <- sqrt(D*(rho2*(1-R22)/(P*(1-P)*J) +
                    (1-rho2)*(1-R12)/(P*(1-P)*J*n)))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.crd2r2",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.crd2r2(mdes=1.391, rho2=.20, n=4, J=20)

cosa.crd2r2 <- function(cn=0, cJ=0, cost=NULL,
                        n=NULL, J=NULL, nJ0=c(10,10),
                        constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                        power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                        rho2, P=.50, RTZ=NULL, k1=-6, k2=6, dist.Z="normal",
                        g2=0, R12=0, R22=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.crd2r2"
  LB <- c(1, g2+2+1)
  df <- quote(J - g2 - 2)
  SSE <- quote(sqrt(D*(rho2*(1-R22)/(P*(1-P)*J) +
                    (1-rho2)*(1-R12)/(P*(1-P)*J*n))))
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

# examples

# cosa given total cost
# cosa.crd2r2(mdes=1.391, cn=1, cJ=10, cost=560, constrain="cost", rho2=.20)

# cosa given per unit costs and power
# cosa.crd2r2(mdes=1.391, cn=1, cJ=10, power=.80, constrain="power", rho2=.20)

# cosa given per unit costs and mdes
# cosa.crd2r2(mdes=1.391, cn=1, cJ=10, power=.80, constrain="mdes", rho2=.20, gm=100)

# cosa (fixed n=10) given total cost
# cosa.crd2r2(mdes=1.391, cn=1, cJ=10, n=10, cost=560, constrain="cost", rho2=.20)

# cosa (fixed n=10) given per unit costs and power
# cosa.crd2r2(mdes=1.391, cn=1, cJ=10, n=10, constrain="power", rho2=.20)

# cosa (fixed n=10) given per unit costs and mdes
# cosa.crd2r2(mdes=1.391, cn=1, cJ=10, n=4, J=20, constrain="mdes", rho2=.20)


