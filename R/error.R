# error handler
.error.handler <- function(x){

  names.x <- names(x)
  if(any(!names.x %in% c("", "n","J","K","L", "g1", "g2", "g3", "g4",
                         "r21","r22","r23","r24", "r2t2", "r2t3", "r2t4",
                         "R12", "R22", "R32", "R42", "RT22", "RT32", "RT42",
                         "rho2", "rho3", "rho4", "omega2", "omega3", "omega4",
                         "q", "q", "p", "P", "alpha", "power", "mdes", "es",
                         "two.tailed", "two.tail"))) {
    stop("Unused arguments", call. = FALSE)
  }

  # exclude default NULL arguments
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
  idx.var <- match(c("r21","r22","r23","r24", "r2t2", "r2t3", "r2t4",
                     "rho2", "rho3", "rho4", "omega2", "omega3", "omega4",
                     "q", "q", "p", "p", "alpha", "power"),  names.x)
  if(any(x[idx.var[!is.na(idx.var)]] < 0) | any(x[idx.var[!is.na(idx.var)]] > 1)){
    stop("Argument out of range when it should be between 0 and 1", call.=FALSE)
  }

  idx.r2 <- match(c("r2", "r21", "r22", "r23", "r24", "r2t2",
                    "r2t3", "r2t4"), names.x)
  if (any(x[idx.r2[!is.na(idx.r2)]] > 0) & any(x[idx.g[!is.na(idx.g)]] == 0)) {
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] > 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] == 0]
    if (any(substr(err.r2, 4, 4) == substr(err.g, 2, 2))) {
      warning("R-squared value for a level may not be greater than zero when the number of covariates at that level [other than blocking variables] is zero",
              call. = FALSE)
    }
  }
  else if (any(x[idx.r2[!is.na(idx.r2)]] == 0) & any(x[idx.g[!is.na(idx.g)]] > 0)) {
    err.r2 <- names.x[idx.r2[!is.na(idx.r2)]][x[idx.r2[!is.na(idx.r2)]] == 0]
    err.g <- names.x[idx.g[!is.na(idx.g)]][x[idx.g[!is.na(idx.g)]] > 0]
    if (any(substr(err.r2, 4, 4) == substr(err.g, 2, 2))) {
      warning("R-squared value for a level may not be zero when the number of covariates at that level [other than blocking variables] is greater than zero ",
              call. = FALSE)
    }
  }

  # warn for negative treatment effect
  if("es" %in% names.x) {
    if(any(x$es < 0)){
      warning("Negative MDES", call.=FALSE)
    }
  }

  # validty check for two-tailed test
  if("two.tailed" %in% names.x){
    if(x$two.tailed != TRUE & x$two.tailed != FALSE){
      stop("Non-logical value for two.tailed", call.=FALSE)
    }
  } else if("two.tail" %in% names.x){
    if(x$two.tail != TRUE & x$two.tail != FALSE){
      stop("Non-logical value for two.tailed", call.=FALSE)
    }
  }


}# .error.handler()
