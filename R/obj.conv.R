# define object conversion functions
optimal.to.mdes <- function(x){
  if(class(x)[1] == "parms" & class(x)[2] == "optimal"){
    idx.par <- match(c("n","J","K","L"), names(x$parms))
    idx.out <- match(c("n","J","K","L"), colnames(x$round.optim))
    for(i in 1:length(idx.par[!is.na(idx.par)])){
      x$parms[idx.par[!is.na(idx.par)][i]] <- x$round.optim[idx.out[!is.na(idx.out)][i]]
    }
    fun.parsed <- scan(text = x$fun, what = "character", sep=".", quiet = TRUE)
    if(length(fun.parsed)==2){
  	  x$fun <- paste0("mdes", ".", fun.parsed[2])
  	}else{
  	  x$fun <- paste0("mdes", ".", fun.parsed[2], ".",  fun.parsed[3])
  	}
    return(do.call(x$fun, x$parms))
  }else{
    stop("x should be an object returned from one of the 'optimal' functions in 'PowerUpR' package.")
  }
}

optimal.to.power <- function(x){
  if(class(x)[1] == "parms" & class(x)[2] == "optimal"){
    idx.par <- match(c("n","J","K","L"), names(x$parms))
    idx.out <- match(c("n","J","K","L"), colnames(x$round.optim))
    for(i in 1:length(idx.par[!is.na(idx.par)])){
      x$parms[idx.par[!is.na(idx.par)][i]] <- x$round.optim[idx.out[!is.na(idx.out)][i]]
    }
    fun.parsed <- scan(text = x$fun, what = "character", sep=".", quiet = TRUE)
    if(length(fun.parsed)==2){
  	  x$fun <- paste0("power", ".", fun.parsed[2])
  	}else{
  	  x$fun <- paste0("power", ".", fun.parsed[2], ".",  fun.parsed[3])
  	}
    return(do.call(x$fun, x$parms))
  }else{
    stop("x should be an object returned from one of the 'optimal' functions in 'PowerUpR' package.")
  }
}

power.to.mdes <- function(x){
  if(class(x)[1] == "parms" & class(x)[2] == "power"){
    x$parms$power <- x$power
    fun.parsed <- scan(text = x$fun, what = "character", sep=".", quiet = TRUE)
    if(length(fun.parsed)==2){
	  x$fun <- paste0("mdes", ".", fun.parsed[2])
	}else{
	  x$fun <- paste0("mdes", ".", fun.parsed[2], ".",  fun.parsed[3])
	}
    return(do.call(x$fun, x$parms))
  }else{
    stop("x should be an object returned from one of the 'power' functions in 'PowerUpR' package.")
  }
}

mdes.to.power <- function(x){
  if(class(x)[1] == "parms" & class(x)[2] == "mdes"){
    x$parms$mdes <- x$mdes[1]
    fun.parsed <- scan(text = x$fun, what = "character", sep=".", quiet = TRUE)
    if(length(fun.parsed)==2){
	    x$fun <- paste0("power", ".", fun.parsed[2])
  	}else{
  	  x$fun <- paste0("power", ".", fun.parsed[2], ".",  fun.parsed[3])
  	}
    return(do.call(x$fun, x$parms, envir=parent.frame()))
  }else{
    stop("x should be an object returned from one of the 'mdes' functions in 'PowerUpR' package.")
  }
}

mdes.to.pctl <- function(x){
  if(class(x)[1] == "parms" & class(x)[2] == "mdes"){
    pctl <- paste0(round(pnorm(x$mdes[1]) * 100, 0), "%")
    pctlL <- paste0(round(pnorm(x$mdes[2]) * 100, 0), "%")
    pctlU <- paste0(round(pnorm(x$mdes[3]) * 100, 0), "%")
    pctlLU <- cbind(pctl, pctlL, pctlU)
    colnames(pctlLU) <- c("pctl", paste0(100 * (1-x$alpha),"% LCL"), paste0(100 * (1-x$alpha),"% LCL"))
    return(pctlLU)
  }else if(is.numeric(x)){
    pctl <- paste0(round(pnorm(x) * 100,0),"%")
    return(pctl)
  }else{
    stop("x should be either an object returned from one of the 'mdes' functions in 'PowerUpR' package or a numeric value.")
  }
}

# defunct functions
mrss.to.mdes <- function(...){
  .Defunct(new = "optimal.to.mdes")
  optimal.to.mdes(...)
}
mrss.to.power <- function(...){
  .Defunct(new = "optimal.to.power")
  optimal.to.power(...)
}
t1t2.error <- function(){
  .Defunct()
}
plot.pars <- function(){
  .Defunct()
}
