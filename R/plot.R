plot.mrss <- plot.power <- plot.mdes <- function(x, ypar = "mdes", xpar = NULL,
                                    xlim = NULL, ylim = NULL,
                                    xlab = NULL, ylab = NULL,
                                    main = NULL, sub = NULL,
                                    locate = FALSE, ...){

  if(any(c("med211", "med221") %in% class(x))) {
    stop("Mediation effects are currently not supported", call. = FALSE)
  }

  if(!ypar %in% c("mdes", "power")){
    stop("Incorrect value for argument 'ypar'", call. = FALSE)
  }

  # design characteristics
  fun.parsed <- scan(text = x$fun, what = "character", sep=".", quiet = TRUE)
  design <- fun.parsed[2]
  nlevels <- substr(design, nchar(design) - 2, nchar(design) - 2)
  if(substr(design, 1, nchar(design) - 3) == "mod") {
    block <- "r"
  } else {
    block <- substr(design, nchar(design) - 1, nchar(design) - 1)
  }

  # default xpar if NULL
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
      stop("Incorrect value for argument 'xpar'", call. = FALSE)
    }
  }

  # object conversions
  capture.output({
    if(ypar == "mdes") {
      if(inherits(x, "mrss"))  x <- mrss.to.mdes(x)
      if(inherits(x, "mdes"))  x <- x
      if(inherits(x, "power"))  x <- power.to.mdes(x)
    } else {
      if(inherits(x, "mrss"))  x <- mrss.to.power(x)
      if(inherits(x, "mdes"))  x <- mdes.to.power(x)
      if(inherits(x, "power"))  x <- x
    }
  })

  # default xlim if NULL
  if(is.null(xlim)) {
    if(substr(design, nchar(design) - 2, nchar(design) - 1) %in% c("22", "33")) {
      xlim <-  c(x$parms$g + 5, 1.5 * x$parms[[xpar]])
    } else {
      xlim <-  c(x$parms$g + 3, 1.5 * x$parms[[xpar]])
    }
  } else {
    if(xlim[1] <= 0 || !is.numeric(xlim) || length(xlim) > 2) {
      stop("Incorrect value for argument 'xlim'", call. = FALSE)
    }
  }

  if(xlim[2] - xlim[1] > 20) {
    xseq <- seq(xlim[1], xlim[2], .25)
  } else {
    xseq <- seq(xlim[1], xlim[2], .125)
  }

  # current values
  names.parms <-  names(x$parms)
  idx <- match(xpar, names.parms)
  x0 <- x$parms[[idx]]
  idy <- match(ypar, names(x))
  y0 <- ifelse(ypar == "power", x[[idy]], x[[idy]][1])

  # plot data
  yout <- matrix(NA, nrow = length(xseq), ncol = 3)
  for(i in 1:nrow(yout)){
    x$parms[idx] <- xseq[i]
    capture.output(
      if(ypar == "mdes"){
        yout[i,] <- do.call(x$fun, x$parms)$mdes
      }else if(ypar == "power"){
        yout[i,1] <- do.call(x$fun, x$parms)$power
      }
    )
  }

  # default ylim if NULL
  if(is.null(ylim)) {
    ifelse(ypar == "mdes", ylim <- range(min(yout[,2]), max(yout[,3])), ylim <- c(0,1))
  }

  # labels
  ifelse(!is.null(ylab), ylab,
         ifelse(ypar == "mdes",
                ifelse(substr(design, 1, nchar(design) - 3) == "mod",
                       ylab <- "Minimum Detectable Effect Size Difference",
                       ylab <- "Minimum Detectable Effect Size"),
                ylab <- "Statistical Power"))


  # plot
  plot.new()
  plot.window(xlim = range(xseq),
              ylim = ylim, ...)
  polygon(c(rev(xseq), xseq), c(rev(yout[,3]), yout[,2]), border = NA,
          col = adjustcolor(4, alpha.f = 0.2))
  lines(xseq, yout[,1], col = adjustcolor(4, alpha.f = 0.5), lty = 1, lwd = 2)
  if(ypar == "mdes") {
    lines(xseq, yout[,2], col = adjustcolor(4, alpha.f = 0.2), lty = 1, lwd = 1.5)
    lines(xseq, yout[,3], col = adjustcolor(4, alpha.f = 0.2), lty = 1, lwd = 1.5)
  }
  title(main = main, sub = sub,
        xlab = ifelse(!is.null(xlab), xlab, xpar),
        ylab = ylab)
  axis(1)
  axis(2)
  box()

  # locate parameters for the current design
  if(locate) {
    points(x0, y0, pch=21, bg = adjustcolor(2, alpha.f = 0.5), cex=1.5)
    abline(v = x0, lty = 5, col = adjustcolor(2, alpha.f = 0.5))
  }

  # benchmark values
  abline(h = ifelse(ypar == "mdes", .20, .80), lty = 5, col = adjustcolor(2, alpha.f = 0.5))

}


