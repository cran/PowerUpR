mdes.bcra3f2 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, P=.50, g2=0, R12=0, R22=0,
                        n, J, K, ...){
  df <- K * (J - 2) - g2
  SSE = sqrt(rho2*(1-R22)/(P*(1-P)*J*K) +
               (1-rho2)*(1-R12)/(P*(1-P)*J*K*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.bcra3f2",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}

# example
# mdes.bcra3f2(rho2=.10, n=20, J=44, K=5)


power.bcra3f2 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, P=.50, g2=0, R12=0, R22=0,
                         n, J, K, ...){
  df <- K * (J - 2) - g2
  SSE <- sqrt(rho2*(1-R22)/(P*(1-P)*J*K) +
                    (1-rho2)*(1-R12)/(P*(1-P)*J*K*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.bcra3f2",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.bcra3f2(rho2=.10, n=20, J=44, K=5)

optimal.bcra3f2 <- function(cn=0, cJ=0, cK=0, cost=NULL, n=NULL, J=NULL, K=NULL,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           nJK0=c(10,10,10), ncase=10, gm=2,
                           constrain="power", optimizer="auglag_cobyla",
                           rho2, P=.50, g2=0, R12=0, R22=0){
  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.bcra3f2"
  LB <- c(1, 3, g2+1)
  df <- quote(K * (J - 2) - g2)
  SSE <- quote(sqrt(rho2*(1-R22)/(P*(1-P)*J*K) +
                (1-rho2)*(1-R12)/(P*(1-P)*J*K*n)))
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

# examples

# optimal sample given total cost
# optimal.bcra3f2(cn=1, cJ=10, cK=100, cost=5600, constrain="cost", rho2=.20)

# optimal sample given per unit costs and power
# optimal.bcra3f2(cn=1, cJ=10, cK=100, power=.80, constrain="power", rho2=.20)

# optimal sample given per unit costs and mdes
# optimal.bcra3f2(cn=1, cJ=10, cK=100, power=.80, constrain="mdes",  rho2=.20)

# conditional optimal sample (fixed n=10) given total cost
# optimal.bcra3f2(cn=1, cJ=10, cK=100, n=10, cost=5600, constrain="cost",  rho2=.20)

# conditional optimal sample (fixed J=10) given per unit costs and power
# optimal.bcra3f2(cn=1, cJ=10, cK=100, J=10, constrain="power", rho2=.20)

# conditional optimal sample (fixed n=10) given per unit costs and mdes
# optimal.bcra3f2(cn=1, cJ=10, cK=100, K=3, constrain="mdes", rho2=.20)

