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
  length.list.n <- length(unlist(x[idx.n[!is.na(idx.n)]]))
  length.unlist.n <- length(x[idx.n[!is.na(idx.n)]])
  if(length.list.n == length.unlist.n){
    if(any(x[idx.n[!is.na(idx.n)]] <= 0)){
      stop("Negative sample size", call.=FALSE)
    }
  }

  # validity check for number of covariates
  idx.g <- match(c("g1", "g2", "g3", "g4"),  names.x)
  if(any(x[idx.g[!is.na(idx.g)]] < 0)){
    stop("Negative number of covariates", call.=FALSE)
  }

  # validity check for variance parameters, proportions, and probabilities
  idx.var <- match(c("R12","R22","R32","R42", "RT22", "RT32", "RT42",
                     "rho2", "rho3", "rho4", "omega2", "omega3", "omega4",
                     "Q", "alpha", "power"),  names.x)
  if(any(x[idx.var[!is.na(idx.var)]] < 0) | any(x[idx.var[!is.na(idx.var)]] > 1)){
    stop("Argument out of range when it should be between 0 and 1", call.=FALSE)
  }
  if(any(x$P < .01 | x$P > .99 )){
    stop("P should be between .01 and .99", call.=FALSE)
  }
  idx.r2 <- match(c("R12","R22","R32","R42", "RT22", "RT32", "RT42"),  names.x)
  if(any(x[idx.r2[!is.na(idx.r2)]] > 0) & any(x[idx.g[!is.na(idx.g)]] == 0)){
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] > 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] == 0]
    if(substr(err.r2, 2, 2) == substr(err.g, 2, 2)){
      stop("R-squared value for a level cannot be greater than zero
           when the number of covariates at that level is zero ", call.=FALSE)
    }
  }else if(any(x[idx.r2[!is.na(idx.r2)]] == 0) & any(x[idx.g[!is.na(idx.g)]] > 0 )){
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] == 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] > 0]
    if(any(substr(err.r2, 2, 2) == substr(err.g, 2, 2))){
      stop("R-squared value for a level cannot be zero
           when the number of covariates at that level is greater than zero ", call.=FALSE)
    }
  }

  # warn for negative treatment effect
  if(any(x$mdes < 0)){
    warning("Negative MDES", call.=FALSE)
  }

  # validty check for two-tailed test
  if(x$two.tail != TRUE & x$two.tail != FALSE){
    stop("Non-logical value for two.tail", call.=FALSE)
  }

  # validty check for other cosa related parameters
  if("cost" %in% names.x){
    if(x$gm <= 0){
      stop("gm should be >= 1", call.=FALSE)
    }
    if(x$ncase < 2){
      stop("ncase should be >= 2", call.=FALSE)
    }
  }

}

