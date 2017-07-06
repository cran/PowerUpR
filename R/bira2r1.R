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

optimal.bira2r1 <- function(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           nJ0=c(10,10), ncase=10, gm=2,
                           constrain="power", optimizer="auglag_cobyla",
                           rho2, omega2, P=.50, g2=0, R12=0, RT22=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.bira2r1"
  LB <- c(1, g2+1+1)
  df <- quote(J - g2 - 1)
  SSE <- quote(sqrt(rho2*omega2*(1-RT22)/J +
                (1-rho2)*(1-R12)/(P*(1-P)*J*n)))
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

# examples

# optimal sample given total cost
# optimal.bira2r1(cn=1, cJ=10, cost=560, constrain="cost", rho2=.20, omega2=.50)

# optimal sample given per unit costs and power
# optimal.bira2r1(cn=1, cJ=10, power=.80, constrain="power", rho2=.20, omega2=.50)

# optimal sample given per unit costs and mdes
# optimal.bira2r1(cn=1, cJ=10, power=.80, constrain="mdes", rho2=.20, omega2=.50)

# conditional optimal sample (fixed n=10) given total cost
# optimal.bira2r1(cn=1, cJ=10, n=10, cost=560, constrain="cost", rho2=.20, omega2=.50)

# conditional optimal sample (fixed n=10) given per unit costs and power
# optimal.bira2r1(cn=1, cJ=10, n=10, constrain="power", rho2=.20, omega2=.50)

# conditional optimal sample (fixed n=10) given per unit costs and mdes
# optimal.bira2r1(cn=1, cJ=10, n=10, constrain="mdes", rho2=.20, omega2=.50)
