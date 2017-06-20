# shiny application to find appropriate function
.PowerUpR.Shiny <- function() {
  shiny::runApp(appDir = system.file("PowerUpR.Shiny", package="PowerUpR"))
}

# minimum detectable effect size
.mdes.fun <- function(power, alpha, SSE, df, two.tail, envir = parent.frame()){
  T1 <- ifelse(two.tail == TRUE, abs(qt(alpha/2, df)), abs(qt(alpha, df)))
  T2 <- abs(qt(power, df))
  M <- envir$ncp <- ifelse(power >= 0.5, T1 + T2, T1 - T2)
  mdes <- M * SSE
  LCI <- mdes * (1 - T1/M)
  UCI <- mdes * (1 + T1/M)
  MLU <- cbind(mdes, LCI, UCI)
  colnames(MLU) <- c("mdes", paste0(100 * (1 - alpha), "% LCL"), paste0(100 * (1 - alpha), "% LCL"))
  return(MLU)
}

# statistical power
.power.fun <- function(mdes, alpha, SSE, df, two.tail, envir = parent.frame()){
  lamda <- envir$ncp <- mdes/SSE
  power <- ifelse(two.tail == FALSE,
                  1 - pt(qt(alpha, df, lower.tail = FALSE), df, lamda),
                  1 - pt(qt(alpha/2, df, lower.tail = FALSE), df, lamda) +
                    pt(-qt(alpha/2, df, lower.tail = FALSE), df, lamda))
  return(power)
}

# constrained optimal sample allocation
.optimal.fun <- function(cn = 0, cJ = 0, cK = 0, cL = 0, cost = NULL, n = NULL, J = NULL, K = NULL, L = NULL,
                         nJ0 = c(10,10), nJK0 = c(10,10,10), nJKL0 = c(10,10,10,10), ncase = 10, gm = 2,
                         constrain = "power", optimizer = "auglag_cobyla",
                         rho2, rho3, rho4, omega2, omega3, omega4,
                         P=.50, Q=NULL, g1=0, g2=0, g3=0, g4=0,
                         R12=0, R22=0, R32=0, R42=0, RT22=0, RT32=0, RT42=0, R2T2=0, R3T2=0,
                         ...){

  df <- get("df", parent.frame())
  SSE <- get("SSE", parent.frame())
  two.tail <- get("two.tail", parent.frame())
  alpha <- get("alpha", parent.frame())
  power <- get("power", parent.frame())
  mdes <- get("mdes", parent.frame())
  LB <- get("LB", parent.frame())
  nlevels <- length(LB)

  # equality constraints on cost
  .eq.cost <- function(ss){
    n <- ss[1]
    J <- ss[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ss[3]
      }
    if(nlevels == 4){
      L <- ss[4]
      }
    cost0 <- switch(as.character(nlevels),
                    '2' = {cJ*J + cn*n*J - cost},
                    '3' = {cK*K + cJ*J*K + cn*n*J*K - cost},
                    '4' = {cL*L + cK*K*L + cJ*J*K*L + cn*n*J*K*L - cost})
    return(cost0)
  }

  # equality constraints on power
  .eq.power <- function(ss){
    n <- ss[1]
    J <- ss[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ss[3]
      }
    if(nlevels == 4){
      L <- ss[4]
      }
    power1 <- .power.fun(mdes, alpha, eval(SSE), eval(df), two.tail)
    power0 <- power1-power
    return(power0)
  }

  # equality constraints on mdes
  .eq.mdes <- function(ss){
    n <- ss[1]
    J <- ss[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ss[3]
      }
    if(nlevels == 4){
      L <- ss[4]
      }
    mdes1 <- .mdes.fun(power, alpha, eval(SSE), eval(df), two.tail)[1]
    mdes0 <-mdes1-mdes
    return(mdes0)
  }

  # minimize cost given power or mdes
  .min.cost <- function(ss){
    n <- ss[1]
    J <- ss[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ss[3]
      }
    if(nlevels == 4){
      L <- ss[4]
      }
    cost <- switch(as.character(nlevels),
                   '2' = {cJ*J + cn*n*J},
                   '3' = {cK*K + cJ*J*K + cn*n*J*K},
                   '4' = {cL*L + cK*K*L + cJ*J*K*L + cn*n*J*K*L})
    return(cost)
  }

  # minimize treatment variance given cost
  .min.var <- function(ss){
    n <- ss[1]
    J <- ss[2]
    if(nlevels == 3 | nlevels == 4){
      K <- ss[3]
      }
    if(nlevels == 4){
      L <- ss[4]
      }
    sigmaT2 <- eval(SSE)^2
    return(sigmaT2)
  }

  fn.constr <- switch(constrain,
                      "power" = ".eq.power",
                      "mdes" = ".eq.mdes",
                      "cost" = ".eq.cost"
                      )

  fn.min <- switch(constrain,
                   "power" = ".min.cost",
                   "mdes" = ".min.cost",
                   "cost" = ".min.var"
                   )

  local.optim <- switch(optimizer,
                        "auglag_cobyla" = "COBYLA",
                        "auglag_lbfgs" = "LBFGS",
                        "auglag_mma" = "MMA",
                        "auglag_slsqp" = "SLSQP"
                        )

  # constrain specified sample sizes, define limits
  if(nlevels == 2){
    n0  <- ifelse(!is.null(n), n, nJ0[1])
    nLB <- ifelse(!is.null(n), n, LB[1])
    nUB <- ifelse(!is.null(n), n, Inf)
    J0  <- ifelse(!is.null(J), J, nJ0[2])
    JLB <- ifelse(!is.null(J), J, LB[2])
    JUB <- ifelse(!is.null(J), J, Inf)
    ss0  <- c(n0,J0)
    ssLB <- c(nLB,JLB)
    ssUB <- c(nUB,JUB)
  }else if(nlevels == 3){
    n0  <- ifelse(!is.null(n), n, nJK0[1])
    nLB <- ifelse(!is.null(n), n, LB[1])
    nUB <- ifelse(!is.null(n), n, Inf)
    J0  <- ifelse(!is.null(J), J, nJK0[2])
    JLB <- ifelse(!is.null(J), J, LB[2])
    JUB <- ifelse(!is.null(J), J, Inf)
    K0  <- ifelse(!is.null(K), K, nJK0[3])
    KLB <- ifelse(!is.null(K), K, LB[3])
    KUB <- ifelse(!is.null(K), K, Inf)
    ss0  <- c(n0,J0,K0)
    ssLB <- c(nLB,JLB,KLB)
    ssUB <- c(nUB,JUB,KUB)
  }else if(nlevels == 4){
    n0  <- ifelse(!is.null(n), n, nJKL0[1])
    nLB <- ifelse(!is.null(n), n, LB[1])
    nUB <- ifelse(!is.null(n), n, Inf)
    J0  <- ifelse(!is.null(J), J, nJKL0[2])
    JLB <- ifelse(!is.null(J), J, LB[2])
    JUB <- ifelse(!is.null(J), J, Inf)
    K0  <- ifelse(!is.null(K), K, nJKL0[3])
    KLB <- ifelse(!is.null(K), K, LB[3])
    KUB <- ifelse(!is.null(K), K, Inf)
    L0  <- ifelse(!is.null(L), L, nJKL0[4])
    LLB <- ifelse(!is.null(L), L, LB[4])
    LUB <- ifelse(!is.null(L), L, Inf)
    ss0  <- c(n0,J0,K0,L0)
    ssLB <- c(nLB,JLB,KLB,LLB)
    ssUB <- c(nUB,JUB,KUB,LUB)
  }

  # optimize
  nlopt.ss <- nloptr::auglag(x0 = ss0, fn = fn.min, heq = fn.constr,
                             localsolver = local.optim, localtol = 1e-8,
                             lower = ssLB, upper = ssUB)
  if(nlopt.ss$convergence<0 | all(nlopt.ss$par == ss0)){
    warning("Solution may not be feasible. Change default settings")
    }

  # round solution
  ss1 <- round(nlopt.ss$par)
  est.cost <- .min.cost(ss1)
  est.mdes <- .eq.mdes(ss1) + mdes
  est.power <- .eq.power(ss1) + power

  round.optim <- cbind(t(ss1), est.cost, est.mdes, est.power)
  colnames(round.optim) <-  c(c("n", "J", "K", "L")[1:nlevels], "cost", "mdes", "power")


  # grid to approximate an integer solution
  lower.grid <- ss1 - rev(1:nlevels)*gm
  lower.grid[lower.grid < LB] <- LB[LB > lower.grid]
  upper.grid <- ss1 + rev(1:nlevels)*gm
  lu.grid <- cbind(lower.grid, upper.grid)
  ss.list <- list()
  for(i in 1:nlevels){
    var <- paste("var", i, sep = "")
    ss.list[[var]] <- lu.grid[i,1]:lu.grid[i,2]
  }
  gridss <- as.matrix(expand.grid(ss.list))

  # difference from specified power (default, .80), mdes (default, .25), or cost
  integer.optim <- matrix(NA,nrow(gridss), nlevels+3)
  integer.optim[,1:nlevels] <- gridss
  for(i in 1:nrow(gridss)){
    integer.optim[i,nlevels+1] <- .min.cost(gridss[i,])
    integer.optim[i,nlevels+2] <- .eq.mdes(gridss[i,])
    integer.optim[i,nlevels+3] <- .eq.power(gridss[i,])
  }

  integer.optim <- round(integer.optim, digits = 3)

  # output result
  if(constrain == "power"){
    idx <- order(abs(integer.optim[,nlevels+3]), integer.optim[,nlevels+1],
                 decreasing = FALSE)[1:ncase]
  }else if(constrain == "mdes"){
    idx <- order(abs(integer.optim[,nlevels+2]), integer.optim[,nlevels+1],
                 decreasing = FALSE)[1:ncase]
  }else if(constrain == "cost"){
    idx <- order(abs(integer.optim[,nlevels+1]-cost), -integer.optim[,nlevels+3],
                 decreasing = FALSE)[1:ncase]
  }

  integer.optim <- integer.optim[idx,]
  integer.optim[,nlevels+2] <- integer.optim[,nlevels+2] + mdes
  integer.optim[,nlevels+3] <- integer.optim[,nlevels+3] + power
  colnames(integer.optim) <- c(c("n", "J", "K", "L")[1:nlevels], "cost", "mdes", "power")

  nloptr <- list(pars = nlopt.ss$par,
                 obj.value = nlopt.ss$value,
                 iter = nlopt.ss$iter,
                 global.solver = nlopt.ss$global_solver,
                 local.solver = nlopt.ss$local_solver,
                 convergence = nlopt.ss$convergence,
                 message = nlopt.ss$message)
  optim.out <- list(fun = get("fun", parent.frame()),
                    parms = get("parms", parent.frame()),
                    nloptr = nloptr,
                    round.optim = round(round.optim, digits = 3),
                    integer.optim = round(integer.optim, digits = 3)
  )

  return(invisible(optim.out))
}

# type I and type II error plot
.t1t2 <- function(ncp, df, alpha, two.tail){

  # t-critical, power, beta
  talpha <- ifelse(two.tail == FALSE,
                   qt(alpha, df, lower.tail = FALSE),
                   qt(alpha/2, df, lower.tail = FALSE)
                   )
  power <- ifelse(two.tail == FALSE,
                  1 - pt(talpha, df, ncp),
                  1 - pt(talpha, df, ncp) +
                      pt(-talpha, df, ncp)
                  )
  beta = 1 - power

  # define functions
  funt0 <- function(x){
    dt(x, df = df, ncp = 0)
    } # central t
  funt1 <- function(x){
    dt(x, df = df, ncp = ncp)
    } #non-central t

  # plot central t distribution
  plot(funt0, xlim = c(-3,10), ylim = c(0, 0.5),
       yaxs = "i", xaxs = "i", bty = "l",
       main = expression(paste("Type I ", (alpha), " and "," Type II ", (beta), " Error Rates")),
       sub = paste("Type I = ", round(alpha ,digits = 3), ",",
                 "Type II = ", round(beta, digits = 3), ",",
                 "ncp =  ", round(ncp, digits = 3)),
       xlab = "", ylab = "")

  par(new = TRUE)

  # plot non-central t distribution
  plot(funt1, xlim = c(-3,10), ylim = c(0, 0.5),
       yaxs = "i", xaxs = "i", yaxt = "n", xaxt = "n",
       bty = "l", xlab = "", ylab = "")
    legend("topright",
         c(expression(alpha),
           expression(beta),
           ifelse(two.tail == TRUE,
                  expression(t[alpha/2]),
                  expression(t[alpha])
                  ),
           expression(ncp)),
         x.intersp = 0.9, y.intersp = 0.9,
         pch = c(19, 19, NA, NA), lty = c(NA, NA, 2, 2), cex = 1,
         col = c(adjustcolor(2, alpha.f = 0.3), adjustcolor(4, alpha.f = 0.3), 2, 4),
         bty = 1)

    # axes labels and subtitle
    title(ylab = "P(t)", line = 2)
    title(xlab = "t", line = 2)

    # draw vertical lines
    abline(v = 0, lty = 2, col = 4) # mean of central t in dashed blue line
    abline(v = ncp, lty = 2, col = 4) # mean of non-central t in dashed blue line
    abline(v = talpha, lty = 2, col = 2) # t-critical in dashed red line

    # shaded area in red for alpha
    xalpha0 <- seq(from = talpha, to = 10, by = .001)
    yalpha0 <- rep(NA, length(xalpha0))
    for(i in 1:length(xalpha0)){
      yalpha0[i] <- funt0(xalpha0[i])
      }

    xalpha <- c(xalpha0, rev(xalpha0))
    yalpha <- c(yalpha0, rep(0, length(yalpha0)))
    polygon(x = xalpha, y = yalpha, col = adjustcolor(2, alpha.f = 0.3))

    # shaded area in light blue for beta
    xbeta0 <- seq(from = -3,to = talpha, by = .001)
    ybeta0 <- rep(NA, length(xbeta0))
    for(i in 1:length(xbeta0)){
      ybeta0[i] <- funt1(xbeta0[i])
      }
    xbeta <- c(xbeta0, rev(xbeta0))
    ybeta <- c(ybeta0, rep(0, length(ybeta0)))
    polygon(x = xbeta, y = ybeta, col = adjustcolor(4, alpha.f = 0.3))
}


# Type I  and Type II error plot wrapper
.t1t2.error <- function(x){

  if(class(x)[2] == "optimal"){
    x <- optimal.to.power(x)
  }

  # plot type I and type II error rates
  return(.t1t2(ncp = x$ncp,
               df = x$parms$df,
               alpha = x$parms$alpha,
               two.tail = x$parms$two.tail))
}

# comparison plots
.ggplot.parms <- function(x, out.parm="power", in.parm, in.seq, by.parm, by.seq, color=TRUE, ...){
  if(class(x)[1] == "parms" & class(x)[2] %in% c("mdes", "power", "optimal")){
    if(out.parm == in.parm){
      stop(paste("'in.parm' and 'out.parm' arguments cannot be same"))
    }
    if(in.parm %in% names(x$parms) == 0){
      stop(paste(in.parm, "is not a valid value for 'in.parm' argument"))
    }
    if(by.parm %in% names(x$parms) == 0){
      stop(paste(by.parm, "is not a valid parameter for 'by.parm' argument"))
    }

    fun.parsed <- scan(text=x$fun, what="character", sep=".", quiet = TRUE)
    out <- data.frame(cbind(0, expand.grid(in.seq, by.seq)))
    colnames(out) <- c("out.parm", "in.parm", "by.parm")

    invisible(capture.output(
      if(out.parm == "power"){
        x <- switch(fun.parsed[1],
                    "mdes" = mdes.to.power(x),
                    "optimal" = optimal.to.power(x),
                    x
        )
        idx.in <- match(in.parm, names(x$parms))
        idx.by <- match(by.parm, names(x$parms))
        for(i in 1:nrow(out)){
          x$parms[c(idx.in, idx.by)] <- out[i, c(2, 3)]
          out[i,1] <- do.call(x$fun, x$parms)$power
        }
      }else if(out.parm == "mdes"){
        x <- switch(fun.parsed[1],
                    "power" = power.to.mdes(x),
                    "optimal" = optimal.to.mdes(x),
                    x
        )
        idx.in <- match(in.parm, names(x$parms))
        idx.by <- match(by.parm, names(x$parms))
        for(i in 1:nrow(out)){
          x$parms[c(idx.in, idx.by)] <- out[i, c(2, 3)]
          out[i,1] <- do.call(x$fun, x$parms)$mdes[1]
        }
      }else{stop(paste(out.parm,"is not a valid parameter for 'out.parm' argument"))}
    ))

    if(color == TRUE){
      # plot sequences
      ggplot2::ggplot(out, aes(in.parm, out.parm, group = as.factor(by.parm), color = as.factor(by.parm))) +
        geom_line() +
        geom_point() +
        xlab(in.parm) + ylab(out.parm) +
        scale_color_discrete(name = by.parm,
                             breaks = as.factor(by.seq),
                             labels = as.factor(by.seq))
    }else{
      ggplot2::ggplot(out, aes(in.parm, out.parm, group = as.factor(by.parm), shape = as.factor(by.parm))) +
        geom_line() +
        geom_point() +
        xlab(in.parm) + ylab(out.parm) +
        scale_shape_discrete(name = by.parm,
                             breaks = as.factor(by.seq),
                             labels = as.factor(by.seq))
    }
  }else{
    stop("x should an object returned from one of the functions in 'PowerUpR' package")
  }
}

# error handler
.error.handler <- function(x){

  # exclude default NULL arguments
  names.x <- names(x)
  idx.notnull <- match(names(lapply(x, is.null)[!lapply(x, is.null) == TRUE]), names.x)
  parms.notnull <- x[idx.notnull]
  # redefine the check list
  names.x <- names(parms.notnull)
  x <-  parms.notnull

  # validity check for sample sizes
  idx.n <- match(c("n","J","K","L"),  names.x)
  if(any(x[idx.n[!is.na(idx.n)]] <= 0)){
    stop(
      simpleError(
        paste("Argument", sQuote(names.x[idx.n[!is.na(idx.n)]][x[idx.n[!is.na(idx.n)]] <= 0]), "cannot be negative")
      )
    )
  }
  if(any(x[idx.n[!is.na(idx.n)]] == 1)){
    warning(
      simpleWarning(
        paste("Argument", sQuote(names.x[idx.n[!is.na(idx.n)]][x[idx.n[!is.na(idx.n)]] == 1]), "is 1")
      )
    )
  }

  # validity check for number of covariates
  idx.g <- match(c("g1", "g2", "g3", "g4"),  names.x)
  if(any(x[idx.g[!is.na(idx.g)]] < 0)){
    stop(
      simpleError(
        paste("Argument", sQuote(names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] < 0]), "cannot be negative")
      )
    )
  }

  # validity check for variance parameters, proportions, and probabilities
  idx.var <- match(c("R12","R22","R32","R42", "RT22", "RT32", "RT42",
                     "rho2", "rho3", "rho4", "omega2", "omega3", "omega4", "P", "Q", "alpha", "power"),  names.x)
  if(any(x[idx.var[!is.na(idx.var)]] < 0) | any(x[idx.var[!is.na(idx.var)]] > 1)){
    stop(
      simpleError(
        paste("Argument", sQuote(names.x[idx.var[!is.na(idx.var)]][x[idx.var[!is.na(idx.var)]] < 0 | x[idx.var[!is.na(idx.var)]] > 1]), "should be between 0 and 1")
      )
    )
  }

  idx.r2 <- match(c("R12","R22","R32","R42", "RT22", "RT32", "RT42"),  names.x)
  if(any(x[idx.r2[!is.na(idx.r2)]] > 0) & any(x[idx.g[!is.na(idx.g)]] == 0)){
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] > 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] == 0]
    if(substr(err.r2, 2, 2) == substr(err.g, 2, 2)){
      warning(
        simpleWarning(
          paste("Value of argument", sQuote(err.r2), "should correspond to the value of argument", sQuote(err.g), ", results may not be precise")
        )
      )
    }
  }else if(any(x[idx.r2[!is.na(idx.r2)]] == 0) & any(x[idx.g[!is.na(idx.g)]] != 0 )){
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] > 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] == 0]
    if(substr(err.r2, 2, 2) == substr(err.g, 2, 2)){
      stop(
        simpleError(
          paste("Value of argument", sQuote(err.r2), "should correspond to the value of argument", sQuote(err.g))
        )
      )
    }
  }

  # warn for negative treatment effect
  if(any(x$mdes < 0)){
    warning(
      simpleWarning(paste("Argument", sQuote(quote(mdes)), "has negative value"))
    )
  }

  # validty check for two-tailed test
  if(x$two.tail != TRUE & x$two.tail != FALSE){
    stop(
      simpleError(paste("Value for argument", sQuote(quote(two.tail)), "should be either 'TRUE' or 1, 'FALSE'or 0"))
    )
  }

  # validty check for other optimal design related parameters
  if("gm" %in% names.x){
    if(x$gm <= 0){
      stop(
        simpleError(paste("Value for argument", sQuote(quote(gm)), "should be >= 1"))
      )
    }
    if(x$ncase < 2){
      stop(
        simpleError(paste("Value for argument", sQuote(quote(ncase)), "should be >= 2"))
      )
    }
  }

}

