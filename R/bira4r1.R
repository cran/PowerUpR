mdes.bira4r1 <- function(power=.80, alpha=.05, two.tail=TRUE,
                        rho2, rho3, rho4, omega2, omega3, omega4,
                        P=.50, R12=0, RT22=0, RT32=0, RT42=0, g4=0,
                        n, J, K, L, ...){
  df <- L - g4 - 1
  SSE <- sqrt(rho4*omega4*(1-RT42)/L +
                rho3*omega3*(1-RT32)/(K*L) +
                rho2*omega2*(1-RT22)/(J*K*L) +
                (1-rho4-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  mdes <- .mdes.fun(power, alpha, SSE, df, two.tail)
  mdes.out <- list(fun = "mdes.bira4r1",
                   parms = parms[-length(parms)],
                   ncp = ncp, mdes = mdes)
  class(mdes.out) <- c("parms", "mdes")
  print(round(mdes.out$mdes, digits=3))
  return(invisible(mdes.out))
}
# example
# mdes.bira4r1(rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, omega2=.50, n=10, J=4, L=27, K=4)


power.bira4r1 <- function(mdes=.25, alpha=.05, two.tail=TRUE,
                         rho2, rho3, rho4, omega2, omega3, omega4,
                         P=.50, R12=0, RT22=0, RT32=0, RT42=0, g4=0,
                         n, J, K, L, ...){
  df <- L - g4 - 1
  SSE <- sqrt(rho4*omega4*(1-RT42)/L +
                       rho3*omega3*(1-RT32)/(K*L) +
                       rho2*omega2*(1-RT22)/(J*K*L) +
                       (1-rho4-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n))
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  power <- .power.fun(mdes, alpha, SSE, df, two.tail)
  power.out <-  list(fun = "power.bira4r1",
                     parms = parms[-length(parms)],
                     ncp = ncp, power = power)
  class(power.out) <- c("parms", "power")
  print(paste("power = ", round(power.out$power, digits=3)))
  return(invisible(power.out))
}
# example
# power.bira4r1(rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, omega2=.50, n=10, J=4, L=27, K=4)

optimal.bira4r1 <- function(cn=0, cJ=0, cK=0, cL=0, cost=NULL, n=NULL, J=NULL, K=NULL, L=NULL,
                           power=.80, mdes=.25, alpha=.05, two.tail=TRUE,
                           nJKL0=c(10,10,10,10), ncase=10, gm=2,
                           constrain="power", optimizer="auglag_cobyla",
                           rho4, rho3, rho2, omega4, omega3, omega2,
                           P=.50, g4=0, RT42=0, RT32=0, RT22=0, R12=0){

  # design-specific
  parms <- as.list(environment(), all=TRUE)
  .error.handler(parms)
  fun = "optimal.bira4r1"
  LB <- c(1, 1, 1, g4+1+1)
  df <- quote(L - g4 - 1)
  SSE <- quote(sqrt(rho4*omega4*(1-RT42)/L +
                rho3*omega3*(1-RT32)/(K*L) +
                rho2*omega2*(1-RT22)/(J*K*L) +
                (1-rho4-rho3-rho2)*(1-R12)/(P*(1-P)*J*K*L*n)))
  # optimal-specific
  optim.out <- do.call(".optimal.fun", parms)
  class(optim.out) <- c("parms", "optimal")
  print(round(optim.out$round.optim, digits=3))
  return(invisible(optim.out))
}

# examples

# optimal sample given total cost
# optimal.bira4r1(cn=1, cJ=10, cK=100, cL=1000, cost=75600, constrain="cost",
#               rho4=.10, rho3=.15, rho2=.20, omega4=.50, omega3=.50, omega2=.50)

# optimal sample given per unit costs and power
# optimal.bira4r1(cn=1, cJ=10, cK=100, cL=1000, power=.80, constrain="power",
#               rho4=.10, rho3=.15, rho2=.20, omega4=.50, omega3=.50, omega2=.50)

# optimal sample given per unit costs and mdes
# optimal.bira4r1(cn=1, cJ=10, cK=100, cL=1000, power=.80, constrain="mdes",
#               rho4=.10, rho3=.15, rho2=.20, omega4=.50, omega3=.50, omega2=.50)

# conditional optimal sample (fixed n=10, and J=2) given total cost
# optimal.bira4r1(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, cost=75600, constrain="cost",
#               rho4=.10, rho3=.15, rho2=.20, omega4=.50, omega3=.50, omega2=.50)

# conditional optimal sample (fixed n=10, and J=2) given per unit costs and power
# optimal.bira4r1(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, power=.80, constrain="power",
#              rho4=.10, rho3=.15, rho2=.20, omega4=.50, omega3=.50, omega2=.50)

# conditional optimal sample (fixed n=10, and J=2) given per unit costs and mdes
# optimal.bira4r1(cn=1, cJ=10, cK=100, cL=1000, n=10, J=2, power=.80, constrain="mdes",
#               rho4=.10, rho3=.15, rho2=.20, omega4=.50, omega3=.50, omega2=.50)