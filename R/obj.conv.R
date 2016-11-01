# define object conversion functions
optimal.to.mdes <- function(design){
  idx.par <- match(c("n","J","K","L"), names(design$par))
  idx.out <- match(c("n","J","K","L"), colnames(design$round.optim))
  for(i in 1:length(idx.par[!is.na(idx.par)])){
    design$par[idx.par[!is.na(idx.par)][i]] <- design$round.optim[idx.out[!is.na(idx.out)][i]]
  }
  design$fun <- paste0("mdes", substr(design$fun,8,14))
  return(do.call(design$fun,design$par))
}
optimal.to.power <- function(design){
  idx.par <- match(c("n","J","K","L"), names(design$par))
  idx.out <- match(c("n","J","K","L"), colnames(design$round.optim))
  for(i in 1:length(idx.par[!is.na(idx.par)])){
    design$par[idx.par[!is.na(idx.par)][i]] <- design$round.optim[idx.out[!is.na(idx.out)][i]]
  }
  design$fun <- paste0("power", substr(design$fun,8,14))
  return(do.call(design$fun,design$par))
}
mrss.to.mdes <- function(design){
  design$fun <- paste0("mdes", substr(design$fun,5,11))
  return(do.call(design$fun,design$par))
}
mrss.to.power <- function(design){
  design$fun <- paste0("power", substr(design$fun,5,11))
  return(do.call(design$fun,design$par))
}
