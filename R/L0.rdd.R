# normal (Schochet, 2008, p.14)
# P: proportion of cases in treatment
.RTZ.n <- function(P){
  RTZ <- dnorm(qnorm(1-P)) / sqrt(P*(1-P))
  return(RTZ)
}

# uniform (Schochet, 2008, p.14)
# P: proportion of cases in treatment
.RTZ.u <- function(P){
  RTZ <- sqrt(3*P*(1-P))
  return(RTZ)
}

# truncated normal (Schochet, 2008, p.14)
# P: proportion of cases in treatment
# k1: left truncation point (in standard deviation units from full normal distribution mean)
# k2: right truncation point (in standard deviation units from full normal distribution mean)
.RTZ.tn <- function(k1, k2, P){
  c <- qnorm(P*pnorm(k1) + (1-P)*pnorm(k2) )
  sigmaZ2 <- 1 - ((k2*dnorm(k2) - k1*dnorm(k1)) / (pnorm(k2) - pnorm(k1))) -
    ((dnorm(k2) - dnorm(k1))/(pnorm(k2) - pnorm(k1)))^2
  RTZ <- (P/sqrt(sigmaZ2*P*(1-P))) * ((dnorm(k2) - dnorm(k1)) / (pnorm(k2) - pnorm(k1)) -
                                        (dnorm(k2) - dnorm(c)) / (pnorm(k2) - pnorm(c)))
  return(RTZ)
}

# design effect
.D.fun <- function(P, k1, k2, dist.Z, RTZ){
  .DE <- function(RTZ){
    D <- 1/(1-RTZ^2)
    return(D)
  }
  if(is.null(RTZ)){
    if(dist.Z=="normal"){
      D <- .DE(.RTZ.tn(k1,k2,P))
    }else if(dist.Z=="uniform"){
      if(k1!=-6 | k2!=6){warning("k1 and/or k2 will be ignored.")}
      D <- .DE(.RTZ.u(P))
    }
  }else if(!is.null(RTZ)){
    D <- .DE(RTZ)
    warning("Make sure RTZ is consistent with P as both are related! See Schochet (2008, p.14).
            With RTZ=0 results are identical to random assignment case.")
  }
  return(D)
  }
