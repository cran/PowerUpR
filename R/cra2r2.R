mdes.cra2r2 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                        rho2, p=.50, g2=0, r21=0, r22=0,
                        n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  df <- J - g2 - 2
  SSE <- sqrt(rho2*(1-r22)/(p*(1-p)*J) +
               (1-rho2)*(1-r21)/(p*(1-p)*J*n))

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "main", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
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

mdesd.mod221 <- mdes.mod221 <- function(power=.80, alpha=.05, two.tailed=TRUE,
                        rho2, omegam2, g1=0, r21=0, r2m2=0,
                        p=.50, q=NULL, n, J){

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(omegam2 == 0 || r2m2 == 1) {
    df <- n*J - J - g1 - 2
    cat("\nNon-randomly varying",
        ifelse(is.null(q), "continuous moderator", "binary moderator"))
  } else {
    df <- J - g1 - 2
    cat("\nRandomly varying",
        ifelse(is.null(q), "continuous moderator", "binary moderator"))
  }

  SSE <- ifelse(is.null(q),
                sqrt(rho2*omegam2*(1-r2m2)/(p*(1-p)*J) +
                       (1-rho2)*(1-r21)/(p*(1-p)*J*n)), # continuous mod
                sqrt(rho2*omegam2*(1-r2m2)/(p*(1-p)*J) +
                       (1-rho2)*(1-r21)/(p*(1-p)*q*(1-q)*J*n)) # binary mod
  )

  mdes <- .mdes.fun(power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed)
  .summ.mdes(effect = "mod", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
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
# examples
# mdes.mod221(omegam2=.2, r2m2=1, rho2=.20, n=4, J=20, q=.3)
# mdes.mod221(rho2=.2, omegam2=0, r2m2=.2, n=4, J=20)
# mdes.mod221(rho2=.2, omegam2=.2, r2m2=.2, q=.3, n=4, J=20)

mdesd.mod222 <- mdes.mod222 <- function(power=.80, alpha=.05, two.tailed=TRUE,
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
  .summ.mdes(effect = "mod", power = power, alpha = alpha, sse = SSE, df = df, two.tailed = two.tailed, mdes = mdes)
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
# examples
# mdes.mod222(rho2=.20, n=4, J=20)
# mdes.mod222(rho2=.20, n=4, J=20, q=.5)

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
        ifelse(is.null(q), "continuous moderator", "binary moderator"))
  } else {
    df <- J - g1 - 2
    cat("\nRandomly varying",
        ifelse(is.null(q), "continuous moderator", "binary moderator"))
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
# examples
# power.mod221(rho2=.2, omegam2=.2, r2m2=.2, n=4, J=20)
# power.mod221(rho2=.2, omegam2=.2, r2m2=.2, q=.3, n=4, J=20)


power.mod222 <- function(es=.25, alpha=.05, two.tailed=TRUE,
                              rho2, g2=0, r21=0, r22=0,
                              p=.50, q=NULL, n, J) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  cat(ifelse(is.null(q),
             "\nContinuous moderator",
             "\nBinary moderator"))

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
# examples
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
    c(round(powera, 3), NA, NA, NA),
    c(round(powerb1, 3), NA, NA, NA),
    c(round(powerB, 3), NA, NA, NA),
    c(NA, round(power.sobel.ab1, 3), round(power.joint.ab1, 3), round(power.mc.ab1, 3)),
    c(NA, round(power.sobel.ab2, 3), round(power.joint.ab2, 3), round(power.mc.ab2, 3)),
    c(NA, round(power.sobel.aB, 3), round(power.joint.aB, 3), round(power.mc.aB, 3))
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
                    power = round(power, 3))
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
    c(round(powera, 3), NA, NA, NA),
    c(round(powerb, 3), NA, NA, NA),
    c(NA, round(power.sobel.ab, 3), round(power.joint.ab, 3), round(power.mc.ab, 3))
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
        ifelse(is.null(q), "continuous moderator", "binary moderator"),
        "\nJ =", J)
  } else {
    cat("\nRandomly varying",
        ifelse(is.null(q), "continuous moderator", "binary moderator"),
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

