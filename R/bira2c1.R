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

optimal.bira2c1 <- function(cn=0, cJ=0, cost=NULL, n=NULL, J=NULL,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           nJ0=c(10,10), ncase=10, gm=2,
                           constrain="power", optimizer="auglag_cobyla",
                           P=.50, g1=0, R12=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.bira2c1"
  LB <- c(2, g1+1+1)
  df <- quote(J * (n - 1) - g1 - 1)
  SSE <- quote(sqrt((1-R12)/(P*(1-P)*J*n)))
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

# examples

# optimal sample given total cost
# optimal.bira2c1(cn=1, cJ=10, cost=560, constrain="cost")

# optimal sample given per unit costs and power
# optimal.bira2c1(cn=1, cJ=10, constrain="power", gm=300)

# optimal sample given per unit costs and mdes
# optimal.bira2c1(cn=1, cJ=10, constrain="mdes", gm=300)

# conditional optimal sample (fixed n=10) given total cost
# optimal.bira2c1(cn=1, cJ=10, n=10, cost=560, constrain="cost")

# conditional optimal sample (fixed n=10) given per unit costs and power
# optimal.bira2c1(cn=1, cJ=10, n=10, constrain="power")

# conditional optimal sample (fixed n=10) given per unit costs and mdes
# optimal.bira2c1(cn=1, cJ=10, n=10, constrain="mdes")