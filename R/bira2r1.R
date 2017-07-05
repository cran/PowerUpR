mdes.bira2r1 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2,  omega2, P=.50, g2=0, R12=0, RT22=0,
                        n, J, ...){
  df <- J - g2 - 1
  SSE <- sqrt(rho2*omega2*(1-RT22)/J +
               (1-rho2)*(1-R12)/(P*(1-P)*J*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.bira2r1",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.bira2r1(rho2=.35, omega2=.10, n=83, J=480)

power.bira2r1 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2,  omega2, g2=0, P=.50, R12=0, RT22=0,
                         n, J, ...){
  df <- J - g2 - 1
  SSE <- sqrt(rho2*omega2*(1-RT22)/J +
                    (1-rho2)*(1-R12)/(P*(1-P)*J*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.bira2r1",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.bira2r1(rho2=.35, omega2=.10, n=83, J=480)

cosa.bira2r1 <- function(cn=0, cJ=0, cost=NULL,
                           n=NULL, J=NULL, P=NULL,
                           nJ0=c(10,10), P0=.50,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                           rho2, omega2, g2=0, R12=0, RT22=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.bira2r1"
  LB <- c(1, g2+1+1)
  df <- quote(J - g2 - 1)
  SSE <- quote(sqrt(rho2*omega2*(1-RT22)/J +
                (1-rho2)*(1-R12)/(P*(1-P)*J*n)))
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

optimal.bira2r1 <- function(...){
  .Deprecated(new="cosa.bira2r1")
  cosa.bira2r1(...)
}

# examples

# cosa given total cost
# cosa.bira2r1(cn=1, cJ=10, cost=560, constrain="cost", rho2=.20, omega2=.50)

# cosa given per unit costs and power
# cosa.bira2r1(cn=1, cJ=10, power=.80, constrain="power", rho2=.20, omega2=.50)

# cosa given per unit costs and mdes
# cosa.bira2r1(cn=1, cJ=10, power=.80, constrain="mdes", rho2=.20, omega2=.50)

# cosa (fixed n=10) given total cost
# cosa.bira2r1(cn=1, cJ=10, n=10, cost=560, constrain="cost", rho2=.20, omega2=.50)

# cosa (fixed n=10) given per unit costs and power
# cosa.bira2r1(cn=1, cJ=10, n=10, constrain="power", rho2=.20, omega2=.50)

# cosa (fixed n=10) given per unit costs and mdes
# cosa.bira2r1(cn=1, cJ=10, n=10, constrain="mdes", rho2=.20, omega2=.50)
