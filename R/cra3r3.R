mdes.cra3r3 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                             rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0,
                             n, J, K){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- K - g3 - 2
  SSE <- sqrt(rho3*(1-r23)/(p*(1-p)*K) +
                rho2*(1-r22)/(p*(1-p)*J*K) +
                (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*K*n))

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "main", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.cra3r3",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3,
                                g3=g3, r21=r21, r22=r22, r23=r23,
                                p=p, n=n, J=J, K=K),
                   df = df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("main", "mdes")
  return(invisible(mdes.out))
}

# example
# mdes.cra3r3(rho3=.06, rho2=.17, n=15, J=3, K=60)

mdesd.mod331 <- mdes.mod331 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                        rho2, rho3, omegam2=0, omegam3=0,
                        g1=0, r21=0, r2m2=0, r2m3=0,
                        p=.50, q=NULL, n, J, K){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(omegam2 == 0 || r2m2 == 1) {
    df <- n*J*K - J*K - K - g1 - 2
    if(omegam3 != 0 || r2m3 != 1) {
      omegam3 <- 0
      r2m3 <- 1
      message("Arguments 'omegam3' and/or 'r2m3' are ignored")
    }
    cat("\nModerator effect: Non-randomly varying (level 2 and 3) \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  } else if(omegam3 == 0 || r2m3 == 1) {
    df <- J*K - K - g1 - 2
    cat("\nModerator effect: Non-randomly varying (level 3) \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  } else {
    df <- K - g1 - 1
    cat("\nModerator effect: Randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  }

  SSE <- ifelse(is.null(q),
                sqrt(rho3*omegam3*(1-r2m3)/(p*(1-p)*K) +
                       rho2*omegam2*(1-r2m2)/(p*(1-p)*K*J) +
                       (1-rho2-rho3)*(1-r21)/(p*(1-p)*K*J*n)), # continuous mod
                sqrt(rho3*omegam3*(1-r2m3)/(p*(1-p)*K) +
                       rho2*omegam2*(1-r2m2)/(p*(1-p)*K*J) +
                       (1-rho2-rho3)*(1-r21)/(p*(1-p)*q*(1-q)*K*J*n)) # binary mod
  )

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "mod", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)

  mdes.out <- list(fun = "mdes.mod331",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3, omegam2=omegam2, omegam3=omegam3,
                                r21=r21, r2m2=r2m2, r2m3=r2m3,
                                g1=g1, p=p, q=q, n=n, J=J, K=K),
                   df = df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("mod", "mdes")
  return(invisible(mdes.out))
}

# examples
# mdes.mod331(rho3=.05, rho2=.12, omegam2=.08, omegam3=.07, p=.4, q=.7, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)
# mdes.mod331(rho3=.05, rho2=.12, omegam2=.08, omegam3=.07, p=.4, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)
# mdes.mod331(rho3=.05, rho2=.12, omegam2=0.1, omegam3=0, p=.4, q=.7, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)
# mdes.mod331(rho3=.05, rho2=.12, omegam2=0, omegam3=0, p=.4, g1=1, r21=.20, r2m2=1, r2m3=1,  n=20, J=4, K=60)

mdesd.mod332 <- mdes.mod332 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                              rho2, rho3, omegam3, g2=0, r21=0, r22=0, r2m3=0,
                              p=.50, q=NULL, n, J, K){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(omegam3 == 0 || r2m3 == 1) {
    df <- J*K - K - g2 - 2
    cat("\nModerator effect: Non-randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  } else {
    df <- K - g2 - 1
    cat("\nModerator effect: Randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  }

  SSE <- ifelse(is.null(q),
                sqrt(rho3*omegam3*(1-r2m3)/(p*(1-p)*K) +
                       rho2*(1-r22)/(p*(1-p)*K*J) +
                       (1-rho2-rho3)*(1-r21)/(p*(1-p)*K*J*n)), # continuous mod
                sqrt(rho3*omegam3*(1-r2m3)/(p*(1-p)*K) +
                       rho2*(1-r22)/(p*(1-p)*q*(1-q)*K*J) +
                       (1-rho2-rho3)*(1-r21)/(p*(1-p)*q*(1-q)*K*J*n)) # binary mod
  )

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "mod", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)

  mdes.out <- list(fun = "mdes.mod332",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3, omegam3=omegam3,
                                r21=r21, r22=r22, r2m3=r2m3,
                                g2=g2, p=p, q=q, n=n, J=J, K=K),
                   df = df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("mod332", "mdes")
  return(invisible(mdes.out))
}

# examples
# mdes.mod332(rho3=.1, rho2=.1, omegam3=.05, q=.5, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4, K=60)
# mdes.mod332(rho3=.1, rho2=.1, omegam3=.05, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4, K=60)
# mdes.mod332(rho3=.1, rho2=.1, omegam3=0, q=.5, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4, K=60)
# mdes.mod332(rho3=.1, rho2=.1, omegam3=0, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4, K=60)

mdesd.mod333 <- mdes.mod333 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                             rho2, rho3,  g3=0, r21=0, r22=0, r23=0,
                             p=.50, q=NULL, n, J, K){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  cat(ifelse(is.null(q), "\nModerator type: Continuous\n", "\nModerator type: Binary\n"))

  df <- K-g3-4
  SSE <- ifelse(is.null(q),
                sqrt(rho3*(1-r23)/(p*(1-p)*(K-g3-4)) +
                       rho2*(1-r22)/(p*(1-p)*J*(K-g3-4)) +
                       (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*(K-g3-4)*n)), # continuous mod
                sqrt(rho3*(1-r23)/(p*(1-p)*q*(1-q)*(K-g3-4)) +
                       rho2*(1-r22)/(p*(1-p)*q*(1-q)*J*(K-g3-4)) +
                       (1-rho3-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*J*(K-g3-4)*n)) # binary mod
  )

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "mod", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)

  mdes.out <- list(fun = "mdes.mod333",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3, r21=r21, r22=r22, r23=r23,
                                g3=g3, p=p, q=q, n=n, J=J, K=K),
                   df = df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("mod333", "mdes")
  return(invisible(mdes.out))
}

# examples
# mdes.mod333(rho3=.1, rho2=.1, q=.5, g3=1, r21=.3, r22=.4, r23=.5, n=20, J=4, K=60)
# mdes.mod333(rho3=.1, rho2=.1, g3=1, r21=.3, r22=.4, r23=.5, n=20, J=4, K=60)

power.cra3r3 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                         rho2, rho3, g3=0, r21=0, r22=0, r23=0,
                         p=.50, n, J, K){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- K - g3 - 2
  SSE <- sqrt(rho3*(1-r23)/(p*(1-p)*K) +
                rho2*(1-r22)/(p*(1-p)*J*K) +
                (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*K*n))

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)
  power.out <-  list(fun = "power.cra3r3",
                     parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                  rho2=rho2, rho3=rho3,
                                  g3=g3, r21=r21, r22=r22, r23=r23,
                                  p=p, n=n, J=J, K=K),
                     df = df,
                     ncp = es/SSE,
                     power = power)
  class(power.out) <- c("main", "power")
  return(invisible(power.out))
}

# example
# power.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3, K=60)

power.mod331 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                               rho2, rho3, omegam2, omegam3,
                               g1=0, r21=0, r2m2=0, r2m3=0,
                               p=.50, q=NULL, n, J, K){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(omegam2 == 0 || r2m2 == 1) {
    df <- n*J*K - J*K - K - g1 - 2
    if(omegam3 != 0 || r2m3 != 1) {
      omegam3 <- 0
      r2m3 <- 1
      message("Arguments 'omegam3' and/or 'r2m3' are ignored")
    }
    cat("\nModerator effect: Non-randomly varying (level 2 and 3) \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  } else if(omegam3 == 0 || r2m3 == 1) {
    df <- J*K - K - g1 - 2
    cat("\nModerator effect: Non-randomly varying (level 3) \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  } else {
    df <- K - g1 - 1
    cat("\nModerator effect: Randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  }

  SSE <- ifelse(is.null(q),
                sqrt(rho3*omegam3*(1-r2m3)/(p*(1-p)*K) +
                       rho2*omegam2*(1-r2m2)/(p*(1-p)*K*J) +
                       (1-rho2-rho3)*(1-r21)/(p*(1-p)*K*J*n)), # continuous mod
                sqrt(rho3*omegam3*(1-r2m3)/(p*(1-p)*K) +
                       rho2*omegam2*(1-r2m2)/(p*(1-p)*K*J) +
                       (1-rho2-rho3)*(1-r21)/(p*(1-p)*q*(1-q)*K*J*n)) # binary mod
  )

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)

  power.out <- list(fun = "power.mod331",
                   parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3, omegam2=omegam2, omegam3=omegam3,
                                r21=r21, r2m2=r2m2, r2m3=r2m3,
                                g1=g1, p=p, q=q, n=n, J=J, K=K),
                   df = df,
                   ncp = es/SSE,
                   power = power)
  class(power.out) <- c("mod331", "power")
  return(invisible(power.out))
}

# examples
# power.mod331(es=.16, rho3=.05, rho2=.12, omegam2=.08, omegam3=.07, p=.4, q=.7, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)
# power.mod331(es=.09, rho3=.05, rho2=.12, omegam2=.08, omegam3=.07, p=.4, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)
# power.mod331(es=.16, rho3=.05, rho2=.12, omegam2=0, omegam3=.07, p=.4, q=.7, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)
# power.mod331(es=.09, rho3=.05, rho2=.12, omegam2=0, omegam3=.07, p=.4, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)
# power.mod331(es=.16, rho3=.05, rho2=.12, omegam2=.08, omegam3=0, p=.4, q=.7, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)
# power.mod331(es=.09, rho3=.05, rho2=.12, omegam2=.08, omegam3=0, p=.4, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4, K=60)


power.mod332 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                               rho2, rho3, omegam3, g2=0, r21=0, r22=0, r2m3=0,
                               p=.50, q=NULL, n, J, K){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(omegam3 == 0 || r2m3 == 1) {
    df <- J*K - K - g2 - 2
    cat("\nModerator effect: Non-randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  } else {
    df <- K - g2 - 1
    cat("\nModerator effect: Randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"))
  }

  SSE <- ifelse(is.null(q),
                sqrt(rho3*omegam3*(1-r2m3)/(p*(1-p)*K) +
                       rho2*(1-r22)/(p*(1-p)*K*J) +
                       (1-rho2-rho3)*(1-r21)/(p*(1-p)*K*J*n)), # continuous mod
                sqrt(rho3*omegam3*(1-r2m3)/(p*(1-p)*K) +
                       rho2*(1-r22)/(p*(1-p)*q*(1-q)*K*J) +
                       (1-rho2-rho3)*(1-r21)/(p*(1-p)*q*(1-q)*K*J*n)) # binary mod
  )

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)

  power.out <- list(fun = "power.mod332",
                   parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3, omegam3=omegam3,
                                g2=g2, r21=r21, r22=r22, r2m3=r2m3,
                                p=p, q=q, n=n, J=J, K=K),
                   df = df,
                   ncp = es/SSE,
                   power = power)
  class(power.out) <- c("mod332", "power")
  return(invisible(power.out))
}

# examples
# power.mod332(es=.22, rho3=.1, rho2=.1, omegam3=.05, q=.5, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4, K=60)
# power.mod332(es=.11, rho3=.1, rho2=.1, omegam3=.05, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4, K=60)
# power.mod332(es=.22, rho3=.1, rho2=.1, omegam3=0, q=.5, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4, K=60)
# power.mod332(es=.11, rho3=.1, rho2=.1, omegam3=0, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4, K=60)

power.mod333 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                              rho2, rho3, g3=0, r21=0, r22=0, r23=0,
                              p=.50, q=NULL, n, J, K){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  cat(ifelse(is.null(q), "\nModerator type: Continuous\n", "\nModerator type: Binary\n"))

  df <- K-g3-4
  SSE <- ifelse(is.null(q),
                sqrt(rho3*(1-r23)/(p*(1-p)*(K-g3-4)) +
                       rho2*(1-r22)/(p*(1-p)*J*(K-g3-4)) +
                       (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*(K-g3-4)*n)), # continuous mod
                sqrt(rho3*(1-r23)/(p*(1-p)*q*(1-q)*(K-g3-4)) +
                       rho2*(1-r22)/(p*(1-p)*q*(1-q)*J*(K-g3-4)) +
                       (1-rho3-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*J*(K-g3-4)*n)) # binary mod
  )

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)

  power.out <- list(fun = "power.mod333",
                   parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3,
                                g3=g3, r21=r21, r22=r22, r23=r23,
                                p=p, q=q, n=n, J=J, K=K),
                   df = df,
                   ncp = es/SSE,
                   power = power)
  class(power.out) <- c("mod333", "power")
  return(invisible(power.out))
}

# examples
# power.mod333(es=.30, rho3=.1, rho2=.1, q=.5, g3=1, r21=.3, r22=.4, r23=.5, n=20, J=4, K=60)
# power.mod333(es=.15, rho3=.1, rho2=.1, g3=1, r21=.3, r22=.4, r23=.5, n=20, J=4, K=60)

power.med321 <- function(esa, esB, two.tailed = TRUE, alpha = .05,
                         mc = FALSE, nsims = 1000, ndraws = 1000,
                         rhom3, rho2, rho3, r2m2, r2m3, r21, r22, r23,
                         p = .50, n, J, K) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  dfa <- K - 5
  dfB <- K - 6
  df <- rbind(dfa, dfB)
  colnames(df) <- "df"
  rownames(df) <- c("a", "B")

  sea321 <- .se.a321(rhom3 = rhom3, r2m2 = r2m2, r2m3 = r2m3, p = p, J = J, K = K)
  seB321 <- .se.B321(rhom3 = rhom3, rho2 = rho2, rho3 = rho3, r2m2 = r2m2, r2m3 = r2m3,
                     r21 = r21, r22 = r22, r23 = r23, p = p, n = n, J = J, K = K)
  ncpa <- esa/sea321
  ncpB <- esB/seB321
  ncp <- rbind(ncpa, ncpB)
  colnames(ncp) <- "ncp"
  rownames(ncp) <- c("a", "B")

  powera <- .power.fun(es = esa, alpha = alpha, sse = sea321, two.tailed = two.tailed, df = dfa)
  powerB <- .power.fun(es = esB, alpha = alpha, sse = seB321, two.tailed = two.tailed, df = dfB)
  power.sobel.aB <- .power.sobel(x = esa, y = esB, sex = sea321, sey = seB321, alpha = alpha, two.tailed = two.tailed)
  power.joint.aB <- .power.jt(x = esa, y = esB, sex = sea321, sey = seB321, alpha = alpha, dfx = dfa, dfy = dfB, two.tailed = two.tailed)
  power.mc.aB <- ifelse(mc, .power.mc(nsims = nsims, ndraws = ndraws, x = esa, y =esB, sex = sea321, sey = seB321, alpha = alpha, two.tailed = two.tailed), NA)

  power <- rbind(
    c(round(powera, 3), NA, NA, NA),
    c(round(powerB, 3), NA, NA, NA),
    c(NA, round(power.sobel.aB, 3), round(power.joint.aB, 3), round(power.mc.aB, 3))
  )
  colnames(power) <- c("t", "sobel", "joint", "mc")
  rownames(power) <- c("a", "B", "aB")

  power.out <- list(fun = "power.med321",
                    parms = list(esa = esa, esB = esB,
                                 two.tailed = two.tailed, alpha = alpha,
                                 mc = mc, nsims = nsims, ndraws = ndraws,
                                 rhom3 = rhom3, rho2 = rho2, rho3 = rho3, r2m2 = r2m2, r2m3 = r2m3,
                                 r21 = r21, r22 = r22, r23 = r23, p = p, n = n, J = J, K = K),
                    df = df,
                    ncp = ncp,
                    power = round(power, 3))
  cat("Statistical power: \n")
  cat("------------------------------------ \n")
  print(power)
  cat("------------------------------------ \n")
  cat("Degrees of freedom for path a:", dfa,
      "\nDegrees of freedom for path B:", dfB,
      "\nStandardized standard error for path a:", round(sea321, 3),
      "\nStandardized standard error for path B:", round(seB321, 3),
      "\nType I error rate:", alpha,
      "\nTwo-tailed test:", two.tailed, "\n")
  class(power.out) <- c("power", "med321")
  return(invisible(power.out))
}
# examples
# power.med321(esa= .49, esB = .30, rhom3 = 0.26, rho2 = .15, rho3 = .20,
#              r2m2 = .07, r2m3 = .17, r21 = .02, r22 = .41, r23 = .38,
#              p = .50, n = 20, J = 4, K = 30)
#
# power.med321(esa= .51, esB = .30, rhom3 = 0.27, rho2 = .15, rho3 = .19,
#              r2m2 = .07, r2m3 = .16, r21 = .02, r22 = .41, r23 = .38,
#              p = .50, n = 20, J = 4, K = 60)

mrss.cra3r3 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                        n, J, K0=10, tol=.10,
                        rho2, rho3, p=.50, g3=0, r21=0, r22=0, r23=0){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <- K0-g3-2
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    K1 <- (M/es)^2 * (rho3*(1-r23)/(p*(1-p)) +
                          rho2*(1-r22)/(p*(1-p)*J) +
                          (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*n))
    K0 <- (K1+K0)/2
    i <- i+1
  }
  K <- ifelse(df>0,round(K0),NA)

  mrss.out <-  list(fun = "mrss.cra3r3",
                    parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                 n=n, J=J, K0=K0, tol=tol,
                                 rho2=rho2, rho3=rho3,
                                 p=p, r21=r21, r22=r22, r23=r23, g3=g3),
                    df = df,
                    ncp = M,
                    K = K)
  class(mrss.out) <- c("main", "mrss")
  cat("K =", K, "\n")
  return(invisible(mrss.out))
}
# example
# mrss.cra3r3(rho3=.06, rho2=.17, n=15, J=3)

mrss.mod331 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                              rho2, rho3, omegam2, omegam3,
                              g1=0, r21=0, r2m2=0, r2m3=0,
                              p=.50, q=NULL, n, J, K0=10, tol=.10){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    if(omegam2 == 0 || r2m2 == 1) {
      df <- n*J*K0 - J*K0 - K0 - g1 - 2
      if(omegam3 != 0 || r2m3 != 1) {
        omegam3 <- 0
        r2m3 <- 1
      }
    } else if(omegam3 == 0 || r2m3 == 1) {
      df <- J*K0 - K0 - g1 - 2
    } else {
      df <- K0 - g1 - 1
    }
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    K1 <- ifelse(is.null(q),
                 (M/es)^2 * (rho3*omegam3*(1-r2m3)/(p*(1-p)) +
                               rho2*omegam2*(1-r2m2)/(p*(1-p)*J) +
                               (1-rho2-rho3)*(1-r21)/(p*(1-p)*J*n)), # continuous mod
                 (M/es)^2 * (rho3*omegam3*(1-r2m3)/(p*(1-p)) +
                               rho2*omegam2*(1-r2m2)/(p*(1-p)*J) +
                               (1-rho2-rho3)*(1-r21)/(p*(1-p)*q*(1-q)*J*n))  # binary mod
    )
    K0 <- (K1+K0)/2
    i <- i+1
  }
  K <- ifelse(df>0,round(K0),NA)

  mrss.out <- list(fun = "mrss.mod331",
                   parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3, omegam2=omegam2, omegam3=omegam3,
                                r21=r21, r2m2=r2m2, r2m3=r2m3,
                                g1=g1, p=p, q=q, n=n, J=J, K0=K0, tol=tol),
                   df = df,
                   ncp = M,
                   K = K)

  if(omegam2 == 0 || r2m2 == 1) {
    if(omegam3 != 0 || r2m3 != 1) {
      message("Arguments 'omegam3' and/or 'r2m3' are ignored")
    }
    cat("\nModerator effect: Non-randomly varying (level 2 and 3) \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"),
        "\nK =", K)
  } else if(omegam3 == 0 || r2m3 == 1) {
    cat("\nModerator effect: Non-randomly varying (level 3) \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"),
        "\nK =", K)
  } else {
    cat("\nModerator effect: Randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"),
        "\nK =", K)
  }


  class(mrss.out) <- c("mod331", "mrss")
  return(invisible(mrss.out))
}

# examples
# mrss.mod331(es=.16, rho3=.05, rho2=.12, omegam2=.08, omegam3=.07, p=.4, q=.7, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4)
# mrss.mod331(es=.09, rho3=.05, rho2=.12, omegam2=.08, omegam3=.07, p=.4, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4)
# mrss.mod331(es=.16, rho3=.05, rho2=.12, omegam2=0, omegam3=.07, p=.4, q=.7, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4)
# mrss.mod331(es=.09, rho3=.05, rho2=.12, omegam2=0, omegam3=.07, p=.4, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4)
# mrss.mod331(es=.16, rho3=.05, rho2=.12, omegam2=.08, omegam3=0, p=.4, q=.7, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4)
# mrss.mod331(es=.09, rho3=.05, rho2=.12, omegam2=.08, omegam3=0, p=.4, g1=1, r21=.20, r2m2=0, r2m3=0,  n=20, J=4)



mrss.mod332 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                              rho2, rho3, omegam3, g2=0, r21=0, r22=0, r2m3=0,
                              p=.50, q=NULL, n, J, K0=10, tol=.10){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    if(omegam3 == 0 || r2m3 == 1) {
      df <- J*K0 - K0 - g2 - 2
    } else {
      df <- K0 - g2 - 1
    }
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    K1 <- ifelse(is.null(q),
                 (M/es)^2 * (rho3*omegam3*(1-r2m3)/(p*(1-p)) +
                               rho2*(1-r22)/(p*(1-p)*J) +
                               (1-rho2-rho3)*(1-r21)/(p*(1-p)*J*n)), # continuous mod
                 (M/es)^2 * (rho3*omegam3*(1-r2m3)/(p*(1-p)) +
                               rho2*(1-r22)/(p*(1-p)*q*(1-q)*J) +
                               (1-rho2-rho3)*(1-r21)/(p*(1-p)*q*(1-q)*J*n)) # binary mod
    )
    K0 <- (K1+K0)/2
    i <- i+1
  }
  K <- ifelse(df>0,round(K0),NA)

  mrss.out <- list(fun = "mrss.mod332",
                   parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3, omegam3=omegam3,
                                g2=g2, r21=r21, r22=r22, r2m3=r2m3,
                                p=p, q=q, n=n, J=J, K0=K0, tol=tol),
                   df = df,
                   ncp = M,
                   K = K)

  if(omegam3 == 0 || r2m3 == 1) {
    cat("\nModerator effect: Non-randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"),
        "\nK =", K)
  } else {
    cat("\nModerator effect: Randomly varying \nModerator type:",
        ifelse(is.null(q), "Continuous\n", "Binary\n"),
        "\nK =", K)
  }

  class(mrss.out) <- c("mod332", "mrss")
  return(invisible(mrss.out))
}

# examples
# mrss.mod332(es=.22, rho3=.1, rho2=.1, omegam3=.05, q=.5, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4)
# mrss.mod332(es=.11, rho3=.1, rho2=.1, omegam3=.05, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4)
# mrss.mod332(es=.22, rho3=.1, rho2=.1, omegam3=0, q=.2, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4)
# mrss.mod332(es=.11, rho3=.1, rho2=.1, omegam3=0, g2=1, r21=.30, r22=.4, r2m3=0,  n=20, J=4)

mrss.mod333 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                             rho2, rho3, g3=0, r21=0, r22=0, r23=0,
                             p=.50, q=NULL, n, J, K0=10, tol=.10){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <- K0 - g3 - 2
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    K1 <- ifelse(is.null(q),
                 g3 + 4 + (M/es)^2 * (rho3*(1-r23)/(p*(1-p)) +
                                        rho2*(1-r22)/(p*(1-p)*J) +
                                        (1-rho3-rho2)*(1-r21)/(p*(1-p)*J*n)), # continuous mod
                 g3 + 4 + (M/es)^2 * (rho3*(1-r23)/(p*(1-p)*q*(1-q)) +
                                        rho2*(1-r22)/(p*(1-p)*q*(1-q)*J) +
                                        (1-rho3-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*J*n)) # binary mod
    )
    K0 <- (K1+K0)/2
    i <- i+1
  }
  K <- ifelse(df>0,round(K0),NA)

  mrss.out <- list(fun = "mrss.mod333",
                   parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, rho3=rho3,
                                g3=g3, r21=r21, r22=r22, r23=r23,
                                p=p, q=q, n=n, J=J, K0=K0, tol=tol),
                   df = df,
                   ncp = M,
                   K = K)
  cat(ifelse(is.null(q), "\nModerator type: Continuous\n", "\nModerator type: Binary\n"),
      "\nK =", K)
  class(mrss.out) <- c("mod333", "mrss")
  return(invisible(mrss.out))
}

# examples
# mrss.mod333(es=.30, rho3=.1, rho2=.1, q=.5, g3=1, r21=.3, r22=.4, r23=.5, n=20, J=4)
# mrss.mod333(es=.15, rho3=.1, rho2=.1, g3=1, r21=.3, r22=.4, r23=.5, n=20, J=4)

