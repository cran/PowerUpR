plot.power <- plot.mdes  <- plot.mrss <- function(x, ypar = "power", xpar, xseq, ...){

  object <- x
  user.parms <- as.list(match.call())
  if(any(c("pars", "type", "mrss.seq", "mdes.seq", "power.seq",
            "left.right.angle", "up.down.angle", "nlevels")%in% names(user.parms))) {
    stop("Defunct parameters", call. = FALSE)
  }

  if(!ypar %in% c("mdes", "power")){
    stop("'ypar' can only be ne of the'power' or 'mdes'")
  }
  if(!xpar %in% c("n","J","K","L")){
    stop("'xpar' can only be one of the'n', 'J', 'K', or 'L'")
  }

  capture.output({
    if(inherits(object, "mrss")){
      if(ypar == "mdes"){
        object <- mrss.to.mdes(object)
      }else if(ypar == "power"){
        object <- mrss.to.power(object)
      }
    }

    if(inherits(object, "mdes") & ypar == "power"){
      object <- mdes.to.power(object)
    }
    if(inherits(object, "power") & ypar == "mdes"){
      object <- power.to.mdes(object)
    }
  })

  if(min(xseq) <= 0) {
    stop("Sample size for 'xpar' should be greater than 1")
  }

  names.parms <-  names(object$parms)
  # current values
  idx <- match(xpar, names.parms)
  x0 <- object$parms[[idx]]
  idy <- match(ypar, names(object))
  y0 <- ifelse(ypar == "power", object[[idy]], object[[idy]][1])

  # plot
  x <- expand.grid(xseq, NA)
  for(i in 1:nrow(x)){
    object$parms[idx] <- x[i,1]
    capture.output(
      if(ypar == "mdes"){
        x[i,2] <- do.call(object$fun, object$parms)$mdes[1]
      }else if(ypar == "power"){
        x[i,2]<- do.call(object$fun, object$parms)$power
      }
      )
  }


  plot(x,
       type="l",
       xlab=xpar, ylab=ypar,
       xlim=range(xseq),
       ...)
  points(x0, y0, col='black', pch=21, bg = "red", cex=1.5)
  abline(v = x0, h = ifelse(ypar == "power", .80, .25), lty = 3)
}


