# define object conversion functions
optimal.to.mdes <- function(design){
  idx.par <- match(c("n","J","K","L"), names(design$par))
  idx.out <- match(c("n","J","K","L"), colnames(design$round.optim))
  for(i in 1:length(idx.par[!is.na(idx.par)])){
    design$par[idx.par[!is.na(idx.par)][i]] <- design$round.optim[idx.out[!is.na(idx.out)][i]]
  }
  design$fun <- paste0("mdes", ".", sub('.*\\.', '', design$fun))
  return(do.call(design$fun,design$par))
}

optimal.to.power <- function(design){
  idx.par <- match(c("n","J","K","L"), names(design$par))
  idx.out <- match(c("n","J","K","L"), colnames(design$round.optim))
  for(i in 1:length(idx.par[!is.na(idx.par)])){
    design$par[idx.par[!is.na(idx.par)][i]] <- design$round.optim[idx.out[!is.na(idx.out)][i]]
  }
  design$fun <- paste0("power", ".", sub('.*\\.', '', design$fun))
  return(do.call(design$fun,design$par))
}

mrss.to.mdes <- function(design){
  design$fun <- paste0("mdes", ".", sub('.*\\.', '', design$fun))
  return(do.call(design$fun,design$par))
}

mrss.to.power <- function(design){
  design$fun <- paste0("power", ".", sub('.*\\.', '', design$fun))
  return(do.call(design$fun,design$par))
}

power.to.mdes <- function(design){
  design$par$power <- design$power
  design$fun <- paste0("mdes", ".", sub('.*\\.', '', design$fun))
  return(do.call(design$fun,design$par))
}

mdes.to.pctl <- function(design){
  pctl <- paste0(round(pnorm(design$mdes[1])*100,0),"%")
  pctlL <- paste0(round(pnorm(design$mdes[2])*100,0),"%")
  pctlU <- paste0(round(pnorm(design$mdes[3])*100,0),"%")
  pctlLU <- cbind(pctl, pctlL, pctlU)
  colnames(pctlLU) <- c("pctl", "95% LCL", "95% UCL")
  return(pctlLU)
}
