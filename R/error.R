.error.handler <- function(x) {

  names.x <- names(x)
  if(any(!names.x %in% c("", "n", "J", "K", "L", "n0", "J0", "K0", "L0", "g1", "g2", "g3", "g4",
                         "r21","r22","r23","r24", "r2t2", "r2t3", "r2t4", "rho2", "rho3", "rho4", "rhom2",
                         "omega2", "omega3", "omega4", "omegam2", "omegam3", "r2m1", "r2m2", "r2m3",
                         "q", "p", "alpha", "power", "mdes", "es", "escp", "esa", "esb", "esb1", "esB",
                         "esab", "esab1", "esaB", "esa0", "esb0", "esb10", "esB0", "powera", "powerb",
                         "powerb1", "powerB", "maxiter", "abstol", "tol", "two.tailed",
                         "mc", "nsims", "ndraws", "rhoab", "rhoab1", "rhoab2", "rhob1b2"))) {
    stop("Unused arguments", call. = FALSE)
  }

  # exclude NULL arguments and redefine the check list
  idx.notnull <- match(names(lapply(x, is.null)[!lapply(x, is.null) == TRUE]),
                       names.x)
  parms.notnull <- x[idx.notnull]
  names.x <- names(parms.notnull)
  x <- lapply(parms.notnull, eval)


  # validity check for sample sizes
  idx.n <- intersect(c("n","J","K","L"),  names.x)
  length.list.n <- length(unlist(x[idx.n]))
  length.unlist.n <- length(x[idx.n])
  if(length.list.n == length.unlist.n){
    if(any(x[idx.n] <= 0) ||
       any(lapply(x[idx.n], function(x)!is.numeric(x)) == TRUE) ||
       any(lapply(x[idx.n], length) > 1)) {
      stop("Incorrect sample size", call.=FALSE)
    }
  }

  # validity check for number of covariates
  idx.g <- intersect(c("g1", "g2", "g3", "g4"),  names.x)
  if(length(idx.g) > 0) {
    if(!is.numeric(x[[idx.g]]) ||
       length(x[idx.g]) > 1 ||
       any(x[idx.g] < 0)) {
      stop("Incorrect number of covariates", call.=FALSE)
    }
  }

  # validity check for variance parameters, proportions, and probabilities
  idx.var <- intersect(c("r21","r22","r23","r24", "r2t2", "r2t3", "r2t4",
                         "r2m1", "r2m2", "r2m3","rhom2", "omegam2", "omegam3",
                         "rho2", "rho3", "rho4", "omega2", "omega3", "omega4",
                         "q", "q", "p", "p", "alpha", "power", "powera",
                         "powerb", "powerb1", "powerB"),  names.x)
  if(any(lapply(x[idx.var], function(x)!is.numeric(x)) == TRUE) ||
     any(lapply(x[idx.var], length) > 1) ||
     any(x[idx.var] < 0) ||
     any(x[idx.var] > 1)) {
    stop("Incorrect value for [0, 1] bounded arguments", call.=FALSE)
  }

  # validity check for correlations
  idx.rho <- intersect(c( "rhoab", "rhoab1", "rhoab2", "rhob1b2"),  names.x)
  if(length(idx.rho) > 0) {
    if(any(lapply(x[idx.rho], function(x)!is.numeric(x)) == TRUE) ||
       any(lapply(x[idx.rho], length) > 1) ||
       any(unlist(x[idx.rho]) > 1 ||
           unlist(x[idx.rho]) < -1)) {
      stop("Incorrect value for correlations", call.=FALSE)
    }
  }

  # validity check for R-squared value and number of covariate consistency
  idx.r2 <- intersect(c("r21", "r22", "r23", "r24", "r2t2",
                        "r2t3", "r2t4", "r2m1", "r2m2", "r2m3"), names.x)
  if(length(idx.g) != 0 & length(idx.r2) != 0) {
    if (any(x[idx.r2] > 0) & x[idx.g] == 0) {
      x.r2 <- x[idx.r2]
      x.g <- x[idx.g]
      err.r2 <- names(x.r2[x.r2 > 0])
      err.g <- names(x.g[x.g == 0])
      if (any(substr(err.r2, nchar(err.r2), nchar(err.r2))== substr(err.g, 2, 2))){
        warning("R-squared value for a level may not be greater than zero",
                call. = FALSE)
      }
    } else if (any(x[idx.r2] == 0) & x[idx.g] > 0) {
      x.r2 <- x[idx.r2]
      x.g <- x[idx.g]
      err.r2 <- names(x.r2[x.r2 == 0])
      err.g <- names(x.g[x.g > 0])
      if (any(substr(err.r2, nchar(err.r2), nchar(err.r2)) == substr(err.g, 2, 2))) {
        warning("R-squared value for a level may not be zero",
                call. = FALSE)
      }
    }
  }

  # validity check for effect size
  idx.es <- intersect(c("escp", "es", "esa", "esb", "esb1", "esB",
                        "esa0", "esb0", "esb10", "esB0"),  names.x)
  if(any(lapply(x[idx.es], function(x)!is.numeric(x)) == TRUE) ||
     any(lapply(x[idx.es], length) > 1) ||
     any(x[idx.es] < 0)) {
    stop("Incorrect value for effect size", call.=FALSE)
  }
  if(any(x[idx.es] > 5)) {
    stop("Extreme value for effect size", call.=FALSE)
  }

  # validty check for two-tailed test
  if("two.tailed" %in% names.x){
    if(!is.logical(x$two.tailed) | length(x$two.tailed) > 1 ){
      stop("Non-logical value for argument 'two.tailed'", call.=FALSE)
    }
  }

  # validty check for maxiter
  if("maxiter" %in% names.x){
    if(length(x$maxiter) > 1 ||
       !is.numeric(x$maxiter) ||
       x$maxiter < 10 ||
       x$maxiter > 5000){
      stop("Incorrect value for argument 'maxiter'", call.=FALSE)
    }
  }

  # validty check for mc
  if("mc" %in% names.x){
    if(length(x$mc) > 1 ||
       !is.logical(x$mc)){
      stop("Incorrect value for argument 'mc'", call.=FALSE)
    }
  }

  # validty check for nsims
  if("nsims" %in% names.x){
    if(length(x$nsims) > 1 ||
       !is.numeric(x$nsims) ||
       x$nsims < 10 ||
       x$nsims > 1e6){
      stop("Incorrect value for argument 'nsims'", call.=FALSE)
    }
  }

  # validty check for ndraws
  if("ndraws" %in% names.x){
    if(length(x$ndraws) > 1 ||
       !is.numeric(x$ndraws) ||
       x$ndraws < 10 ||
       x$ndraws > 1e6){
      stop("Incorrect value for argument 'ndraws'", call.=FALSE)
    }
  }

}#.error.handler
