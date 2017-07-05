mdes.bira2c1 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        P=.50, g1=0, R12=0, n, J, ...){
  df <- J * (n - 1) - g1 - 1
  SSE <- sqrt((1-R12)/(P*(1-P)*J*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.bira2c1",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.bira2c1(n=55, J=14)

power.bira2c1 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                        P=.50, g1=0, R12=0, n, J, ...){
  df <- J * (n - 1) - g1 - 1
  SSE <- sqrt((1-R12)/(P*(1-P)*J*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.bira2c1",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.bira2c1(n=55, J=14)

cosa.bira2c1 <- function(cn=0, cJ=0, cost=NULL,
                           n=NULL, J=NULL, P=NULL,
                           nJ0=c(10,10), P0=.50,
                           constrain="power", optimizer="auglag_slsqp", gm=2, ncase=10,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           g1=0, R12=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun <- "cosa.bira2c1"
  LB <- c(2, g1+1+1)
  df <- quote(J * (n - 1) - g1 - 1)
  SSE <- quote(sqrt((1-R12)/(P*(1-P)*J*n)))
  # cosa-specific
  optim.out <- do.call(".cosa.fun", parms)
  class(optim.out) <- c("parms", "cosa")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

optimal.bira2c1 <- function(...){
  .Deprecated(new="cosa.bira2c1")
  cosa.bira2c1(...)
}

# examples

# cosa given total cost
# cosa.bira2c1(cn=1, cJ=10, cost=560, constrain="cost")

# cosa given per unit costs and power
# cosa.bira2c1(cn=1, cJ=10, constrain="power", gm=300)

# cosa given per unit costs and mdes
# cosa.bira2c1(cn=1, cJ=10, constrain="mdes", gm=300)

# cosa (fixed n=10) given total cost
# cosa.bira2c1(cn=1, cJ=10, n=10, cost=560, constrain="cost")

# cosa (fixed n=10) given per unit costs and power
# cosa.bira2c1(cn=1, cJ=10, n=10, constrain="power")

# cosa (fixed n=10) given per unit costs and mdes
# cosa.bira2c1(cn=1, cJ=10, n=10, constrain="mdes")
