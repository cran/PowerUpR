plot.power <- plot.mdes <- function(x, ypar = "power", xpar = NULL,
                                                  ylim = NULL, xlim = NULL, ...){

  object <- x
  names.parms <-  names(object$parms)
  if(any(c("med211", "med221") %in% class(object))) {
    stop("Indirect effects are currently not supported", call. = FALSE)
  }

  user.parms <- as.list(match.call())
  if(any(c("pars", "mrss.seq", "mdes.seq", "power.seq",
            "left.right.angle", "up.down.angle", "nlevels")%in% names(user.parms))) {
    stop("Defunct parameters", call. = FALSE)
  }

  if(!ypar %in% c("mdes", "power")){
    stop("incorrect value for argument 'ypar'", call. = FALSE)
  }

  fun.parsed <- scan(text = object$fun, what = "character", sep=".", quiet = TRUE)
  design <- fun.parsed[2]
  nlevels <- substr(design, nchar(design) - 2, nchar(design) - 2)
  block <- substr(design, nchar(design) - 1, nchar(design) - 1)

  if(is.null(xpar)) {
    if(block == "r") {
      xpar <- switch(nlevels,
                     "1" = "n",
                     "2" = "J",
                     "3" = "K",
                     "4" = "L")
    } else {
      xpar <- switch(nlevels,
                     "2" = "n",
                     "3" = "J",
                     "4" = "K")
    }
  } else {
    if(!xpar %in% c("n","J","K","L")[1:nlevels]){
      stop("Incorrect value for argument 'xpar'")
    }
  }



  if(inherits(object, "mdes") & ypar == "power"){
    capture.output({
      object <- mdes.to.power(object)
    })
  } else if(inherits(object, "power") & ypar == "mdes"){
    capture.output({
      object <- power.to.mdes(object)
    })
  } else if(inherits(object, "mrss")){
    stop("Objects returned from MRSS functions are not allowed", call. = FALSE)
  }

  # current values
  idx <- match(xpar, names.parms)
  x0 <- object$parms[[idx]]
  idy <- match(ypar, names(object))
  y0 <- ifelse(ypar == "power", object[[idy]], object[[idy]][1])

  if(is.null(xlim)) {
    ifelse(block == "r",
           xlim <- switch(nlevels,
                          "1" = c(object$parms$g1 + 3, 1.5 * object$parms[[xpar]]),
                          "2" = c(object$parms$g2 + 3, 1.5 * object$parms[[xpar]]),
                          "3" = c(object$parms$g3 + 3, 1.5 * object$parms[[xpar]]),
                          "4" = c(object$parms$g4 + 3, 1.5 * object$parms[[xpar]])),
           xlim <- switch(nlevels,
                          "2" = c(object$parms$g1 + 3, 1.5 * object$parms[[xpar]]),
                          "3" = c(object$parms$g2 + 3, 1.5 * object$parms[[xpar]]),
                          "4" = c(object$parms$g3 + 3, 1.5 * object$parms[[xpar]]))
    )
  } else {
    if(xlim[1] <= 0 || !is.numeric(xlim) || length(xlim) > 2) {
      stop("Incorrect value for argument 'xlim'")
    }
  }

  # plot
  xseq <- seq(xlim[1], xlim[2], .5)
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

  if(is.null(ylim)) {
    ifelse(ypar == "mdes", ylim <- range(x[,2]), ylim <- c(0,1))
  }

  plot(x, xlab = xpar, ylab = ypar,
       xlim = xlim, ylim = ylim,
       ...)
  points(x0, y0, col='black', pch=21, bg = "red", cex = 1.5)
  abline(v = x0, h = ifelse(ypar == "power", .80, 0), lty = 3)
}


