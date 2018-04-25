mdes.cra2r2 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                        rho2, p=.50, g2=0, r21=0, r22=0,
                        n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- J - g2 - 2
  SSE <- sqrt(rho2*(1-r22)/(p*(1-p)*J) +
               (1-rho2)*(1-r21)/(p*(1-p)*J*n))

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.cra2r2",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, g2=g2, r21=r21, r22=r22,
                                p=p, n=n, J=J),
                   df=df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("main", "mdes")
  return(invisible(mdes.out))
}
# example
# mdes.cra2r2(rho2=.20, n=4, J=20, alpha=.01)

mdes.mod221 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                        rho2, omegam2, g1=0, r21=0, r2m2=0,
                        p=.50, q=NULL, n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(omegam2 == 0 || r2m2 == 1) {
    df <- n*J - J - g1 - 2
    cat("\nNon-randomly varying",
        ifelse(is.null(q), "continuous moderator \n", "binary moderator \n"))
  } else {
    df <- J - g1 - 2
    cat("\nRandomly varying",
        ifelse(is.null(q), "continuous moderator \n", "binary moderator \n"))
  }

  SSE <- ifelse(is.null(q),
                sqrt(rho2*omegam2*(1-r2m2)/(p*(1-p)*J) +
                       (1-rho2)*(1-r21)/(p*(1-p)*J*n)), # continuous mod
                sqrt(rho2*omegam2*(1-r2m2)/(p*(1-p)*J) +
                       (1-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*J*n)) # binary mod
  )

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.mod221",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, omegam2=omegam2, g1=g1, r21=r21, r2m2=r2m2,
                                p=p, q=q, n=n, J=J),
                   df=df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("mod221", "mdes")
  return(invisible(mdes.out))
}
# example
# mdes.mod221(omegam2=.2, r2m2=1, rho2=.20, n=4, J=20, q=.3)
# mdes.mod221(rho2=.2, omegam2=0, r2m2=.2, n=4, J=20)
# mdes.mod221(rho2=.2, omegam2=.2, r2m2=.2, q=.3, n=4, J=20)


mdes.mod222 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                             rho2, g2=0, r21=0, r22=0,
                             p=.50, q=NULL, n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  cat(ifelse(is.null(q), "\nContinuous moderator", "\nBinary moderator"))

  df <- J - g2 - 4
  SSE <- ifelse(is.null(q),
                sqrt(rho2*(1-r22)/(p*(1-p)*(J-g2-4)) +
                       (1-rho2)*(1-r21)/(p*(1-p)*(J-g2-4)*n)), # continuous mod
                sqrt(rho2*(1-r22)/(p*(1-p)*q*(1-q)*(J-g2-4)) +
                       (1-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*(J-g2-4)*n)) # binary mod
  )

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
  mdes.out <- list(fun = "mdes.mod222",
                   parms = list(power=power, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, g2=g2, r21=r21, r22=r22,
                                p=p, q=q, n=n, J=J),
                   df=df,
                   ncp = mdes[1]/SSE,
                   mdes = mdes)
  class(mdes.out) <- c("mod222", "mdes")
  return(invisible(mdes.out))
}
# example
# mdes.mod222(rho2=.20, n=4, J=20)
# mdes.mod222(rho2=.20, n=4, J=20, q=.5)

mdes.med211 <- function(esa0 = .35, esB0 = .35, esb10 = .35,
                        maxiter = 100, abstol = 1e-3,
                        powera = .90, powerb1 = .90, powerB = .90,
                        two.tailed = TRUE, alpha = .05,
                        escp, rhom2, rho2, r22 = 0, r21 = 0,
                        r2m1 = 0, r2m2 = 0, p, n, J) {

  # disable correlated multivariate normal draws
  rhoab1 = 0; rhoab2 = 0; rhob1b2 = 0;

  # argument checks
  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  dfa <- J - 4
  dfb1 <- n * J - 3
  dfB <- J - 5
  df <- rbind(dfa, dfb1, dfB)
  colnames(df) <- "df"
  rownames(df) <- c("a", "b1", "B")

  capture.output({
    # minimum detectable effect size for path a
    ncp <- conv <- i <- 0
    while(i <= maxiter & conv == 0){
      sea0 <- .se.a211(esa = esa0, rhom2 = rhom2, r2m1 = r2m1, r2m2 = r2m2, n = n, J = J, p = p)
      esa1 <- .mdes.fun(power = powera, alpha = alpha, sse = sea0, df = dfa, two.tailed = two.tailed)
      sea1 <- .se.a211(esa = esa1[1], rhom2 = rhom2, r2m1 = r2m1, r2m2 = r2m2, n = n, J = J, p = p)
      powera1 <- .power.fun(es = esa1[1], alpha = alpha, sse = sea1, two.tailed = two.tailed, df = dfa)
      if(abs(powera1 - powera) < abstol){conv <- 1}
      esa0 <- (esa0 + esa1[1])/2
      i <- i+1
    }
    mdesa <- as.vector(esa1)
    names(mdesa) <- c("mdesa", paste0(100 * (1 - alpha), "% lcl"), paste0(100 * (1 - alpha), "% ucl"))

    # minimum detectable effect size for path b1
    ncp <- conv <- i <- 0
    while(i <= maxiter & conv == 0){
      seb10 <- .se.b1211(esb1 =  esb10, rho2 = rho2, rhom2 = rhom2, r21 = r21, r2m1 = r2m1, n = n, J = J)
      esb11 <- .mdes.fun(power = powerb1, alpha = alpha, sse = seb10, df= dfb1, two.tailed = two.tailed)
      seb11 <- .se.b1211(esb1 =  esb11[1], rho2 = rho2, rhom2 = rhom2, r21 = r21, r2m1 = r2m1, n = n, J = J)
      powerb11 <- .power.fun(es = esb11[1], alpha = alpha, sse = seb11 , two.tailed = two.tailed, df = dfb1)
      if(abs(powerb11 - powerb1) < abstol){conv <- 1}
      esb10 <- (esb10 + esb11[1])/2
      i <- i+1
    }
    mdesb1 <- as.vector(esb11)
    names(mdesb1) <- c("mdesb1", paste0(100 * (1 - alpha), "% lcl"), paste0(100 * (1 - alpha), "% ucl"))

    # minimum detectable effect size for path B = b1 + b2
    ncp <- conv <- i <- 0
    while(i <= maxiter & conv == 0){
      seB0 <- .se.B211(esa = esa1[1], esB =  esB0, esb1 =  esb11[1], escp = escp, rho2 = rho2, rhom2 = rhom2, r22 = r22,
                       r21 = r21, r2m2 = r2m2, r2m1 = r2m1, n = n, J = J, p = p)
      esB1 <- .mdes.fun(power = powerB, alpha = alpha, sse = seB0, df= dfB, two.tailed = two.tailed)
      seB1 <- .se.B211(esa = esa1[1], esB = esB1[1], esb1 = esb11[1], escp = escp, rho2 = rho2, rhom2 = rhom2, r22 = r22,
                       r21 = r21, r2m2 = r2m2, r2m1 = r2m1, n = n, J = J, p = p)
      powerB1 <- .power.fun(es = esB1[1], alpha = alpha, sse = seB1 , two.tailed = two.tailed, df = dfB)
      if(abs(powerB1 - powerB) < abstol){conv <- 1}
      esB0 <- (esB0 + esB1[1])/2
      i <- i+1
    }
    mdesB <- as.vector(esB1)
    names(mdesB) <- c("mdesB", paste0(100 * (1 - alpha), "% lcl"), paste0(100 * (1 - alpha), "% ucl"))
  })

  seb21 <- sqrt(seB1^2 + seb11^2)

  ncpa <- mdesa[1]/sea1
  ncpb1 <- mdesb1[1]/seb11
  ncpB <- mdesB[1]/seB1
  ncp <- rbind(ncpa, ncpb1, ncpB)
  colnames(ncp) <- "ncp"
  rownames(ncp) <- c("a", "b1", "B")

  indir.effect.cl <- .mdes3.mc.cl(
    x = mdesa[1], y = mdesb1[1], z = mdesB[1] - mdesb1[1],
    sex = sea1, sey = seb11, sez = seb21,
    rhoxy = rhoab1, rhoxz = rhoab2, rhoyz = rhob1b2,
    alpha = alpha, two.tailed = two.tailed
  )

  mdes <- rbind(mdesa, mdesb1, mdesB, indir.effect.cl)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - alpha), "% lcl"), paste0(100 * (1 - alpha), "% ucl"))
  rownames(mdes) <- c("a", "b1", "B", "ab1", "ab2", "aB")

  mdes.out <- list(fun = "mdes.med211",
                   parms = list(esa0=esa0, esB0=esB0, esb10=esb10,
                                # rhoab1=rhoab1, rhoab2=rhoab2, rhob1b2=rhob1b2,
                                maxiter=maxiter, abstol=abstol,
                                powera=powera, powerb1=powerb1, powerB=powerB,
                                two.tailed=two.tailed, alpha=alpha,
                                escp=escp, rhom2=rhom2, rho2=rho2, r22=r22, r21=r21,
                                r2m1=r2m1, r2m2=r2m2, p=p, n=n, J=J),
                   df = df,
                   ncp = ncp,
                   mdes = mdes)
  cat("Minimum detectable effect size: \n")
  cat("------------------------------------ \n")
  print(mdes)
  cat("------------------------------------ \n")
  cat("Degrees of freedom for path a:", dfa,
      "\nDegrees of freedom for path b1:", dfb1,
      "\nDegrees of freedom for path B:", dfB,
      "\nStandardized standard error for path a:", round(sea1, 3),
      "\nStandardized standard error for path b1:", round(seb11, 3),
      "\nStandardized standard error for path B:", round(seB1, 3),
      "\nType I error rate:", alpha,
      "\nTwo-tailed test:", two.tailed)
  class(mdes.out) <- c("med211", "mdes")
  return(invisible(mdes.out))
}

# example
# dm <- mdes.med211(escp=.1, powera = .90, powerB = .90, powerb1 = .90, maxiter = 500,
#                        rhom2=.3, rho2=.3, r22=.6, r21=.6, r2m2=.6, r2m1=.6, n=30, J=80, p=.1)
# opt.mdes.med(dm, indir.effect = "ab1")

mdes.med221 <- function(esa0 = .35, esb0 = .35,
                               maxiter = 100, abstol = 1e-3,
                               powera = .90, powerb = .90,
                               two.tailed = TRUE, alpha = .05,
                               escp, rho2, r22, r21, r2m2,
                               p, n, J) {

  # disable correlated multivariate normal draws
  rhoab = 0

  # argument checks
  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  dfa <- dfb <- J - 4
  df <- rbind(dfa, dfb)
  colnames(df) <- "df"
  rownames(df) <- c("a", "b")

  capture.output({
    # minimum detectable effect size for path a
    ncp <- conv <- i <- 0
    while(i <= maxiter & conv == 0){
      sea0 <- .se.a221(esa = esa0, r2m2 = r2m2, p = p, J = J)
      esa1 <- .mdes.fun(power = powera, alpha = alpha, sse = sea0, df = dfa, two.tailed = two.tailed)
      sea1 <- .se.a221(esa = esa1[1], r2m2 = r2m2, p = p, J = J)
      powera1 <- .power.fun(es = esa1[1], alpha = alpha, sse = sea1, two.tailed = two.tailed, df = dfa)
      if(abs(powera1 - powera) < abstol){conv <- 1}
      esa0 <- (esa0 + esa1[1])/2
      i <- i+1
    }
    mdesa <- as.vector(esa1)
    names(mdesa) <- c("mdesa", paste0(100 * (1 - alpha), "% lcl"), paste0(100 * (1 - alpha), "% ucl"))

    # minimum detectable effect size for path b
    ncp <- conv <- i <- 0
    while(i <= maxiter & conv == 0){
      seb0 <- .se.b221(esa = esa1[1], esb = esb0, escp = escp, rho2 = rho2, r22 = r22, r21 = r21, r2m2 = r2m2, p = p, n = n, J = J)
      esb1 <- .mdes.fun(power = powerb, alpha = alpha, sse = seb0, df= dfb, two.tailed = two.tailed)
      seb1 <- .se.b221(esa = esa1[1], esb = esb1[1], escp = escp, rho2 = rho2, r22 = r22, r21 = r21, r2m2 = r2m2, p = p, n = n, J = J)
      powerb1 <- .power.fun(es = esb1[1], alpha = alpha, sse = seb1 , two.tailed = two.tailed, df = dfb)
      if(abs(powerb1 - powerb) < abstol){conv <- 1}
      esb0 <- (esb0 + esb1[1])/2
      i <- i+1
    }
    mdesb <- as.vector(esb1)
    names(mdesb) <- c("mdesb", paste0(100 * (1 - alpha), "% lcl"), paste0(100 * (1 - alpha), "% ucl"))
  })

  ncpa <- mdesa[1]/sea1
  ncpb <- mdesb[1]/seb1
  ncp <- rbind(ncpa, ncpb)
  colnames(ncp) <- "ncp"
  rownames(ncp) <- c("a", "b")

  indir.effect.cl <- .mdes2.mc.cl(
    x = mdesa[1], y = mdesb[1],
    sex = sea1, sey = seb1,
    rhoxy = rhoab,
    alpha = alpha, two.tailed = two.tailed
  )

  mdes <- rbind(mdesa, mdesb, indir.effect.cl)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - alpha), "% lcl"), paste0(100 * (1 - alpha), "% ucl"))
  rownames(mdes) <- c("a", "b", "ab")

  mdes.out <- list(fun = "mdes.med221",
                   parms = list(esa0=esa0, esb0=esb0, #rhoab=rhoab,
                                maxiter=maxiter, abstol=abstol,
                                powera=powera, powerb=powerb,
                                two.tailed=two.tailed, alpha=alpha,
                                escp = escp, rho2=rho2, r22=r22, r21=r21, r2m2=r2m2,
                                p=p, n=n, J=J),
                   df = df,
                   ncp = ncp,
                   mdes = mdes)
  cat("Minimum detectable effect size: \n")
  cat("------------------------------------ \n")
  print(mdes)
  cat("------------------------------------ \n")
  cat("Degrees of freedom for path a:", dfa,
      "\nDegrees of freedom for path b:", dfb,
      "\nStandardized standard error for path a:", round(sea1, 3),
      "\nStandardized standard error for path b:", round(seb1, 3),
      "\nType I error rate:", alpha,
      "\nTwo-tailed test:", two.tailed)
  class(mdes.out) <- c("med221", "mdes")
  return(invisible(mdes.out))
}

# example
# dm <- mdes.med221(escp=.1, powera = .90, powerb = .90, rho2=.15, r22=.52, r21=.40, r2m2=.50, n=100, J=40, p=.5)
# opt.mdes.med(dm, indir.effect = "ab")


power.cra2r2 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                         rho2, g2=0, p=.50, r21=0, r22=0,
                         n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- J - g2 - 2
  SSE <- sqrt(rho2*(1-r22)/(p*(1-p)*J) +
                    (1-rho2)*(1-r21)/(p*(1-p)*J*n))

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)
  power.out <-  list(fun = "power.cra2r2",
                     parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                  rho2=rho2,
                                  p=p, r21=r21, r22=r22, g2=g2,
                                  n=n, J=J),
                     df=df,
                     ncp = es/SSE,
                     power = power)
  class(power.out) <- c("main", "power")
  return(invisible(power.out))
}
# example
# power.cra2r2(rho2=.20, n=4, J=20)


power.mod221 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                               rho2, omegam2, g1=0, r21=0, r2m2=0,
                               p=.50, q=NULL, n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(omegam2 == 0 || r2m2 == 1) {
    df <- n*J - J - g1 - 2
    cat("\nNon-randomly varying",
        ifelse(is.null(q), "continuous moderator \n", "binary moderator \n"))
  } else {
    df <- J - g1 - 2
    cat("\nRandomly varying",
        ifelse(is.null(q), "continuous moderator \n", "binary moderator \n"))
  }

  SSE <- ifelse(is.null(q),
                sqrt(rho2*omegam2*(1-r2m2)/(p*(1-p)*J) +
                       (1-rho2)*(1-r21)/(p*(1-p)*J*n)), # continuous mod
                sqrt(rho2*omegam2*(1-r2m2)/(p*(1-p)*J) +
                       (1-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*J*n)) # binary mod
  )

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)

  power.out <- list(fun = "power.mod221",
                   parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, omegam2=omegam2, g1=g1, r21=r21, r2m2=r2m2,
                                p=p, q=q, n=n, J=J),
                   df=df,
                   ncp = es/SSE,
                   power = power)
  class(power.out) <- c("mod221", "power")
  return(invisible(power.out))
}
# example
# power.mod221(rho2=.2, omegam2=.2, r2m2=.2, n=4, J=20)
# power.mod221(rho2=.2, omegam2=.2, r2m2=.2, q=.3, n=4, J=20)


power.mod222 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                              rho2, g2=0, r21=0, r22=0,
                              p=.50, q=NULL, n, J) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  cat(ifelse(is.null(q),
             "\n Continuous moderator \n",
             "\n Binary moderator \n"))

  df <- J - g2 - 4
  SSE <- ifelse(is.null(q),
                sqrt(rho2*(1-r22)/(p*(1-p)*(J-g2-4)) +
                       (1-rho2)*(1-r21)/(p*(1-p)*(J-g2-4)*n)), # continuous mod
                sqrt(rho2*(1-r22)/(p*(1-p)*q*(1-q)*(J-g2-4)) +
                       (1-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*(J-g2-4)*n)) # binary mod
  )

  power <- .power.fun(es = es, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.power(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, es = es)

  power.out <- list(fun = "power.mod222",
                   parms = list(es=es, alpha=alpha, two.tailed=two.tailed,
                                rho2=rho2, g2=g2, r21=r21, r22=r22,
                                p=p, q=q, n=n, J=J),
                   df=df,
                   ncp = es/SSE,
                   power = power)
  class(power.out) <- c("mod222", "power")
  return(invisible(power.out))
}
# example
# power.mod222(rho2=.20, n=4, J=20)
# power.mod222(rho2=.20, q=.3, n=4, J=20)

power.med211 <- function(esa, esb1, esB, escp, two.tailed = TRUE, alpha = .05,
                         mc = FALSE, nsims = 1000, ndraws = 1000,
                         rhom2, rho2, r21, r22, r2m1, r2m2,
                         p, n, J) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)


  dfa <- J - 4
  dfb1 <- n * J - 3
  dfB <- J - 5
  df <- rbind(dfa, dfb1, dfB)
  colnames(df) <- "df"
  rownames(df) <- c("a", "b1", "B")

  sea211 <- .se.a211(esa = esa, rhom2 = rhom2, r2m1 = r2m1, r2m2 = r2m2, n = n, J = J, p = p)
  seb1211 <- .se.b1211(esb1 =esb1, rho2 = rho2, rhom2 = rhom2, r21 = r21, r2m1 = r2m1, n = n, J = J)
  seB211 <- .se.B211(esa = esa, esB = esB, esb1 = esb1, escp = escp, rho2 = rho2, rhom2 = rhom2, r22 = r22,
                     r21 = r21, r2m2 = r2m2, r2m1 = r2m1, n = n, J = J, p = p)
  ncpa <- esa/sea211
  ncpb1 <- esb1/seb1211
  ncpB <- esB/seB211
  ncp <- rbind(ncpa, ncpb1, ncpB)
  colnames(ncp) <- "ncp"
  rownames(ncp) <- c("a", "b1", "B")

  seb2211 <- sqrt(seB211^2 + seb1211^2)

  powera <- .power.fun(es = esa, alpha = alpha, sse = sea211, two.tailed = two.tailed, df = dfa)
  powerb1 <- .power.fun(es = esb1, alpha = alpha, sse = seb1211, two.tailed = two.tailed, df = dfb1)
  powerB <- .power.fun(es = esB, alpha = alpha, sse = seB211, two.tailed = two.tailed, df = dfB)
  power.sobel.ab1 <- .power.sobel(x = esa, y =esb1, sex = sea211, sey = seb1211, alpha = alpha, two.tailed = two.tailed)
  power.sobel.ab2 <- .power.sobel(x = esa, y =esB-esb1, sex = sea211, sey = seb2211, alpha = alpha, two.tailed = two.tailed)
  power.sobel.aB <- .power.sobel(x = esa, y =esB, sex = sea211, sey = seB211, alpha = alpha, two.tailed = two.tailed)
  power.joint.ab1 <- .power.jt(x = esa, y =esb1, sex = sea211, sey = seb1211, alpha = alpha, dfx = dfa, dfy = dfb1, two.tailed = two.tailed)
  power.joint.ab2 <- .power.jt(x = esa, y =esB-esb1, sex = sea211, sey = seb2211, alpha = alpha, dfx = dfa, dfy = dfB, two.tailed = two.tailed)
  power.joint.aB <- .power.jt(x = esa, y =esB, sex = sea211, sey = seB211, alpha = alpha, dfx = dfa, dfy = dfB, two.tailed = two.tailed)
  power.mc.ab1 <- ifelse(mc, .power.mc(nsims = nsims, ndraws = ndraws, x = esa, y =esb1, sex = sea211, sey = seb1211, alpha = alpha, two.tailed = two.tailed), NA)
  power.mc.ab2 <- ifelse(mc, .power.mc(nsims = nsims, ndraws = ndraws, x = esa, y =esB-esb1, sex = sea211, sey = seb2211, alpha = alpha, two.tailed = two.tailed), NA)
  power.mc.aB <- ifelse(mc, .power.mc(nsims = nsims, ndraws = ndraws, x = esa, y =esB, sex = sea211, sey = seB211, alpha = alpha, two.tailed = two.tailed), NA)

  power <- rbind(
    c(powera, NA, NA, NA),
    c(powerb1, NA, NA, NA),
    c(powerB, NA, NA, NA),
    c(NA, power.sobel.ab1, power.joint.ab1, power.mc.ab1),
    c(NA, power.sobel.ab2, power.joint.ab2, power.mc.ab2),
    c(NA, power.sobel.aB, power.joint.aB, power.mc.aB)
  )
  colnames(power) <- c("t", "sobel", "joint", "mc")
  rownames(power) <- c("a", "b1", "B", "ab1", "ab2", "aB")

  power.out <- list(fun = "power.med211",
                    parms = list(esa=esa, esb1=esb1, esB=esB, escp=escp,
                                 two.tailed=two.tailed, alpha=alpha,
                                 mc=mc, nsims=nsims, ndraws=ndraws,
                                 rho2=rho2, rhom2=rhom2, r21=r21, r22=r22, r2m2=r2m2,
                                 p=p, n=n, J=J),
                    df = df,
                    ncp = ncp,
                    power = power)
  cat("Statistical power: \n")
  cat("------------------------------------ \n")
  print(power)
  cat("------------------------------------ \n")
  cat("Degrees of freedom for path a:", dfa,
      "\nDegrees of freedom for path b1:", dfb1,
      "\nDegrees of freedom for path B:", dfB,
      "\nStandardized standard error for path a:", round(sea211, 3),
      "\nStandardized standard error for path b1:", round(seb1211, 3),
      "\nStandardized standard error for path B:", round(seB211, 3),
      "\nType I error rate:", alpha,
      "\nTwo-tailed test:", two.tailed)
  class(power.out) <- c("power", "med211")
  return(invisible(power.out))
}

# example
# power.med211(esa=.6, esB=.5, esb1=.1, escp=.1, rhom2=.3, rho2=.3, r22=.6, r21=.6, r2m2=.6, r2m1=.6, n=30, J=80, p=.1)

# statistical power
power.med221 <- function(esa, esb, escp, two.tailed = TRUE, alpha = .05,
                         mc = FALSE, nsims = 1000, ndraws = 1000,
                         rho2, r22, r21, r2m2,
                         p = .50, n, J) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  dfa <- dfb <- J - 4
  df <- rbind(dfa, dfb)
  colnames(df) <- "df"
  rownames(df) <- c("a", "b")

  sea221 <- .se.a221(esa = esa, r2m2 = r2m2, p = p, J = J)
  seb221 <- .se.b221(esa = esa, esb = esb, escp = escp, rho2 = rho2, r22 = r22, r21 = r21, r2m2 = r2m2, p = p, n = n, J = J)

  ncpa <- esa/sea221
  ncpb <- esb/seb221
  ncp <- rbind(ncpa, ncpb)
  colnames(ncp) <- "ncp"
  rownames(ncp) <- c("a", "b")

  powera <- .power.fun(es = esa, alpha = alpha, sse = sea221, two.tailed = two.tailed, df = dfa)
  powerb <- .power.fun(es = esb, alpha = alpha, sse = seb221, two.tailed = two.tailed, df = dfb)
  power.sobel.ab <- .power.sobel(x = esa, y =esb, sex = sea221, sey = seb221, alpha = alpha, two.tailed = two.tailed)
  power.joint.ab <- .power.jt(x = esa, y =esb, sex = sea221, sey = seb221, alpha = alpha, dfx = dfa, dfy = dfb, two.tailed = two.tailed)
  power.mc.ab <- ifelse(mc, .power.mc(nsims = nsims, ndraws = ndraws, x = esa, y =esb, sex = sea221, sey = seb221, alpha = alpha, two.tailed = two.tailed), NA)
  power <- rbind(
    c(powera, NA, NA, NA),
    c(powerb, NA, NA, NA),
    c(NA, power.sobel.ab, power.joint.ab, power.mc.ab)
  )
  colnames(power) <- c("t", "sobel", "joint", "mc")
  rownames(power) <- c("a", "b", "ab")

  power.out <- list(fun = "power.med221",
                    parms = list(esa=esa, esb=esb, escp=escp, two.tailed=two.tailed, alpha=alpha,
                                 mc=mc, nsims=nsims, ndraws=ndraws,
                                 rho2=rho2, r22=r22, r21=r21, r2m2=r2m2,
                                 p=p, n=n, J=J),
                    df = df,
                    ncp = ncp,
                    power = power)
  cat("Statistical power: \n")
  cat("------------------------------------ \n")
  print(power)
  cat("------------------------------------ \n")
  cat("Degrees of freedom for path a:", dfa,
      "\nDegrees of freedom for path b:", dfb,
      "\nStandardized standard error for path a:", round(sea221, 3),
      "\nStandardized standard error for path b:", round(seb221, 3),
      "\nType I error rate:", alpha,
      "\nTwo-tailed test:", two.tailed)
  class(power.out) <- c("power", "med221")
  return(invisible(power.out))
}
# example
# power.med221(esa=.3, esb=.1, escp=.1, rho2=.3, r22=.6, r21=.6, r2m2=.6, n=30, J=80, p=.1, mc=TRUE)

opt.mdes.med <- function(x, indir.effect, power.joint = .80) {

  if(inherits(x, "med221")) {
    if(indir.effect != "ab") {
      stop("Invalid indirect effect", call. = FALSE)
    }
  }else if(inherits(x, "med211")) {
    if(!indir.effect %in% c("ab1", "aB")) {
      stop("Invalid indirect effect", call. = FALSE)
    }
  } else {
    stop("'x' should be an object returned from MDES functions for indirect effects", call. = FALSE)
  }

  .min.es <- function(pwr) {
    x$parms$powera <- pwr
    switch(indir.effect,
           "ab" = x$parms$powerb <- power.joint / pwr,
           "ab1" = x$parms$powerb1 <- power.joint / pwr,
           "aB" =  x$parms$powerB <- power.joint / pwr)
    x2 <- do.call(paste("mdes", class(x)[1], sep = "."), x$parms)
    esa <- x2$mdes[rownames(x2$mdes) == "a", colnames(x2$mdes) == "mdes"]
    if(inherits(x, "med221")) {
      esb <- x2$mdes[rownames(x2$mdes) == "b", colnames(x2$mdes) == "mdes"]
    } else if(inherits(x, "med211")) {
      switch(indir.effect,
             "ab1" = esb1 <- x2$mdes[rownames(x2$mdes) == "b1", colnames(x2$mdes) == "mdes"],
             "aB" =  esB <- x2$mdes[rownames(x2$mdes) == "B", colnames(x2$mdes) == "mdes"])
    }
    esjoint <- switch(indir.effect,
                      "ab" = abs(esa * esb),
                      "ab1" = abs(esa * esb1),
                      "aB" = abs(esa * esB))
    return(esjoint)
  }

  capture.output({
    optim.out <- stats::optimize(.min.es, interval = c(power.joint + .01, .99))
  })

  x$parms$powera <- optim.out$minimum
  if(inherits(x, "med221")) {
    x$parms$powerb <- power.joint / x$parms$powera
  } else if(inherits(x, "med211")) {
    switch(indir.effect,
           "ab1" = x$parms$powerb1 <- power.joint / x$parms$powera,
           "aB" =  x$parms$powerB <- power.joint / x$parms$powera )
  }

  switch(indir.effect,
         "ab" = cat("\nOptimal power for path a:", x$parms$powera,
                    "\nOptimal power for path b:", x$parms$powerb, "\n"),
         "ab1" = cat("\nOptimal power for path a:", x$parms$powera,
                     "\nOptimal power for path b1:", x$parms$powerb1, "\n"),
         "aB" = cat("\nOptimal power for path a:", x$parms$powera,
                    "\nOptimal power for path B:", x$parms$powerB, "\n"))

  x3 <- do.call(paste("mdes", class(x)[1], sep = "."), x$parms)

  return(invisible(x3))
}


mrss.cra2r2 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                        n, J0=10, tol=.10,
                        rho2, g2=0, p=.50, r21=0, r22=0) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <- J0-g2-2
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    J1 <- (M/es)^2 * (rho2*(1-r22)/(p*(1-p)) +
                          (1-rho2)*(1-r21)/(p*(1-p)*n))
    if(abs(J1-J0)<tol){conv <- TRUE}
    J0 <- (J1+J0)/2
    i <- i+1
  }
  J <- ifelse(df>0,round(J0),NA)

  mrss.out <-  list(fun = "mrss.cra2r2",
                    parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                 n=n, J0=J0, tol=tol,
                                 rho2=rho2,
                                 p=p, r21=r21, r22=r22, g2=g2),
                    df=df,
                    ncp = M,
                    J = J)
  class(mrss.out) <- c("main", "mrss")
  cat("J = ", J, "\n")
  return(invisible(mrss.out))
}
# example
# mrss.cra2r2(rho2=.20, n=4)


mrss.mod221 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                              n, J0=10, tol=.10, rho2, omegam2, g1=0, r21=0, r2m2=0,
                              p=.50, q=NULL) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    if(omegam2 == 0 || r2m2 == 1) {
      df <- n*J0 - J0 - g1 - 2
    } else {
      df <- J0 - g1 - 2
    }
    if(df <= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    J1 <- ifelse(is.null(q),
                 (M/es)^2 * (rho2*omegam2*(1-r2m2)/(p*(1-p)) +
                   (1-rho2)*(1-r21)/(p*(1-p)*n)), # continuous mod
                 (M/es)^2 * (rho2*omegam2*(1-r2m2)/(p*(1-p)) +
                   (1-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*n)) # binary mod
    )

    if(abs(J1-J0)<tol){conv <- TRUE}
    J0 <- (J1+J0)/2
    i <- i+1
  }
  J <- ifelse(df>0,round(J0),NA)

  mrss.out <-  list(fun = "mrss.mod221",
                    parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                 n=n, J0=J0, tol=tol, rho2=rho2, omegam2=omegam2,
                                 g1=g1, r21=r21, r2m2=r2m2, p=p, q=q),
                    df=df,
                    ncp = M,
                    J = J)
  if(omegam2 == 0 || r2m2 == 1) {
    cat("\nNon-randomly varying",
        ifelse(is.null(q), "continuous moderator \n", "binary moderator"),
        "\nJ =", J)
  } else {
    cat("\nRandomly varying",
        ifelse(is.null(q), "continuous moderator \n", "binary moderator"),
        "\nJ =", J)
  }
  class(mrss.out) <- c("mod221", "mrss")
  return(invisible(mrss.out))
}
# example
# mrss.mod221(rho2=.20, omegam2=0, n=4)


mrss.mod222 <- function(es=.25, power=.80, alpha=.05, two.tailed=TRUE,
                             n, J0=10, tol=.10, rho2, g2=0, r21=0, r22=0,
                             p=.50, q=NULL) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  i <- 0
  conv <- FALSE
  while(i<=100 & conv==FALSE){
    df <-  J0 - g2 - 4
    if(df<= 0 | is.infinite(df)){break}
    T1 <- ifelse(two.tailed==TRUE,abs(qt(alpha/2,df)),abs(qt(alpha,df)))
    T2 <- abs(qt(power,df))
    M <- ifelse(power>=.5,T1+T2,T1-T2)
    J1 <- ifelse(is.null(q),
                 g2 + 4 + (M/es)^2 * (rho2*(1-r22)/(p*(1-p)) +
                   (1-rho2)*(1-r21)/(p*(1-p)*n)), # continuous mod
                 g2 + 4 + (M/es)^2 * (rho2*(1-r22)/(p*(1-p)*q*(1-q)) +
                   (1-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*n)) # binary mod
    )

    if(abs(J1-J0)<tol){conv <- TRUE}
    J0 <- (J1+J0)/2
    i <- i+1
  }
  J <- ifelse(df>0,round(J0),NA)

  mrss.out <-  list(fun = "mrss.mod22",
                    parms = list(es=es, power=power, alpha=alpha, two.tailed=two.tailed,
                                 n=n, J0=J0, tol=tol, rho2=rho2, g2=g2, r21=r21, r22=r22,
                                 p=p, q=q),
                    df=df,
                    ncp = M,
                    J = J)
  class(mrss.out) <- c("mod22", "mrss")
  cat(ifelse(is.null(q), "\nContinuous moderator", "\nBinary moderator"),
      "\nJ =", J)
  return(invisible(mrss.out))
}
# example
# mrss.mod222(rho2=.20, n=4)

