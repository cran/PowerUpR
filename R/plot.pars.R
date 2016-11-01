plot.pars <- function(x, pars=c("power","mdes",NA), type="p",
                        left.right.angle=30, up.down.angle=30, nlevels=10,
                        mdes.seq=NULL, power.seq=NULL, mrss.seq=NULL, ...){
  design <-x

  # change object returned from 'optimal.design' and 'mrss.design'
  if(substr(design$fun,1,7)=="optimal"){
    if(pars[1]=="mdes"){
      design <- optimal.to.mdes(design)
    }else if(pars[1]=="power"){
      design <- optimal.to.power(design)
    }
  }else if(substr(design$fun,1,4)=="mrss"){
    if(pars[1]=="mdes"){
      design <- mrss.to.mdes(design)
    }else if(pars[1]=="power"){
      design <- mrss.to.power(design)
    }
  }

  if(length(pars)==3){ # plot three parameters
    # select top level sample size as a third parameter by default
    if(is.na(pars[3])){
      match.mrss  <- match(c("n","J","K","L"), names(design$par))
      idxmrss <- sort(match.mrss, decreasing=TRUE)[1]
      pars[3] <- names(design$par)[idxmrss]
    }
    # validity check
    if(all(pars%in%c("mdes","power","n","J","K","L"))==FALSE){
      stop("Parameters should be 'mdes','power', and one of the 'n','J','K'or'L'.")
    }
    if(type%in%c("p","c")==FALSE){
      stop("'type' should be 'p' for perspective, or 'c' for contour.")
    }
    if(!is.null(mdes.seq)&!is.null(power.seq)&!is.null(mrss.seq)){
      stop("Two of the sequences should be defined at most.")
    }

    # sequences for plotting
    if(substr(design$fun,1,4)=="mdes" & "mdes"%in%pars){
      if(!is.null(mdes.seq)){
        stop("Sequence for 'mdes' is not allowed.")
      }
      z <- "mdes"
      x <- pars[pars!="mdes"][1]
      y <- pars[pars!="mdes"][2]
      idx <- match(c("n","J","K","L"),pars)
      if(x==pars[idx[!is.na(idx)]]){
        if(is.null(mrss.seq)){
          if(substr(design$fun,1,4)!="mrss"){
            idxx <- match(pars[idx[!is.na(idx)]], names(design$par))
            xseq <- seq(max(3, design$par[[idxx]]-5), design$par[[idxx]]+20, 1)
          }else{
            idxxx <- match(pars[idx[!is.na(idx)]],colnames(design$round.mrss))
            xseq <- seq(max(3, design$round.mrss[idxxx]-5), design$round.mrss[idxxx]+20, 1)
          }
        }else{
          xseq <- mrss.seq
        }
      }else if(x=="power"){
        if(is.null(power.seq)){
          if(!is.null(design$par$power)){
            xseq <- seq(max(.20, design$par$power-.50), min(.99, design$par$power+.20), .01)
          }else{
            xseq <- seq(max(.20, design$power-.50), min(.99, design$power+.20), .01)
          }
        }else{
          xseq <- power.seq
        }
      }
      if(y==pars[idx[!is.na(idx)]]){
        if(is.null(mrss.seq)){
          if(substr(design$fun,1,4)!="mrss"){
            idxx <- match(pars[idx[!is.na(idx)]], names(design$par))
            yseq <- seq(max(3, design$par[[idxx]]-5), design$par[[idxx]]+20, 1)
          }else{
            idxxx <- match(pars[idx[!is.na(idx)]],colnames(design$round.mrss))
            yseq <- seq(max(3, design$round.mrss[idxxx]-5), design$round.mrss[idxxx]+20, 1)
          }
        }else{
          yseq <- mrss.seq
        }
      }else if(y=="power"){
        if(is.null(power.seq)){
          if(!is.null(design$par$power)){
            yseq <- seq(max(.20, design$par$power-.50), min(.99, design$par$power+.20), .01)
          }else{
            yseq <- seq(max(.20, design$power-.50), min(.99, design$power+.20), .01)
          }
        }else{
          yseq <- power.seq
        }
      }
    }
    if(substr(design$fun,1,5)=="power" & "power"%in%pars){
      if(!is.null(power.seq)){
        stop("Sequence for 'power' is not allowed.")
      }
      z <- "power"
      x <- pars[pars!="power"][1]
      y <- pars[pars!="power"][2]
      idx <- match(c("n","J","K","L"),pars)
      if(x==pars[idx[!is.na(idx)]]){
        if(is.null(mrss.seq)){
          if(substr(design$fun,1,4)!="mrss"){
            idxx <- match(pars[idx[!is.na(idx)]], names(design$par))
            xseq <- seq(max(3, design$par[[idxx]]-5), design$par[[idxx]]+20, 1)
          }else{
            idxxx <- match(pars[idx[!is.na(idx)]],colnames(design$round.mrss))
            xseq <- seq(max(3, design$round.mrss[idxxx]-5), design$round.mrss[idxxx]+20, 1)
          }
        }else{
          xseq <- mrss.seq
        }
      }else if(x=="mdes"){
        if(is.null(mdes.seq)){
          if(!is.null(design$par$mdes)){
            xseq <- seq(max(.05, design$par$mdes-.50), design$par$mdes+.20, .01)
          }else{
            xseq <- seq(max(.05, design$mdes[1]-.50), design$mdes[1]+.20, .01)
          }
        }else{
          xseq <- mdes.seq
        }
      }
      if(y==pars[idx[!is.na(idx)]]){
        if(is.null(mrss.seq)){
          if(substr(design$fun,1,4)!="mrss"){
            idxx <- match(pars[idx[!is.na(idx)]], names(design$par))
            yseq <- seq(max(3, design$par[[idxx]]-5), design$par[[idxx]]+20, 1)
          }else{
            idxxx <- match(pars[idx[!is.na(idx)]],colnames(design$round.mrss))
            yseq <- seq(max(3, design$round.mrss[idxxx]-5), design$round.mrss[idxxx]+20, 1)
          }
        }else{
          yseq <- mrss.seq
        }
      }else if(y=="mdes"){
        if(is.null(mdes.seq)){
          if(!is.null(design$par$mdes)){
            yseq <- seq(max(.05, design$par$mdes-.50), design$par$mdes+.20, .01)
          }else{
            yseq <- seq(max(.05, design$mdes[1]-.50), design$mdes[1]+.20, .01)
          }
        }else{
          yseq <- mdes.seq
        }
      }
    }

    # indexes for plotting
    idx <- match(x,names(design$par))
    idy <- match(y,names(design$par))
    x0 <- which(xseq==round(design$par[[idx]], digits=2))
    y0 <- which(yseq==round(design$par[[idy]], digits=2))
    idz <- match(z,names(design))

    # z axis
    zout <- matrix(NA,length(xseq),length(yseq))
    for(i in 1:length(xseq)){
      for(j in 1:length(yseq)){
        par <- design$par
        par[[idx]] <- xseq[i]
        par[[idy]] <- yseq[j]
        zout[i,j] <- do.call(design$fun, par)[[idz]][1]
      }
    }

    # plots
    if(type=="p"){
      plot3D <- persp(xseq, yseq, zout,
                      theta=left.right.angle, phi=up.down.angle, shade=.40,
                      ticktype = "detailed",
                      xlab=x, ylab=y, zlab=z,
                      xlim=range(xseq),
                      ylim=range(yseq),
                      zlim=range(zout),
                      ...)
      points(trans3d(xseq[x0], yseq[y0], zout[x0,y0], plot3D), col='black', pch=21, bg = "red", cex=1.5)
    }else if(type=="c"){
      contour(xseq, yseq, zout, nlevels=nlevels,
              xlab=x, ylab=y, ...)
      cont <- contourLines(xseq, yseq, zout, levels=round(zout[x0,y0],digits=2))
      lines(cont[[1]]$x, cont[[1]]$y, lwd=2)
      points(xseq[x0], yseq[y0], col='black', pch=21, bg = "red", cex=1.5)
    }
  }else if(length(pars)==2){# plot two parameters
    # validity check
    if(all(pars%in%c("mdes","power","n","J","K","L"))==FALSE){
      stop("Parameters should be 'mdes' or 'power', and one of the 'n','J','K'or'L'.")
    }
    if(!is.null(mdes.seq)&!is.null(power.seq)&!is.null(mrss.seq)){
      stop("One of the sequences should be defined at most.")
    }

    # sequences for plotting
    if(substr(design$fun,1,4)=="mdes" & "mdes"%in%pars){
      if(!is.null(mdes.seq)){
        stop("Sequence for 'mdes' is not allowed.")
      }
      y <- "mdes"
      x <- pars[pars!="mdes"][1]
      idx <- match(c("n","J","K","L"),pars)

      if(any(!is.na(idx))){
        if(is.null(mrss.seq)){
          if(substr(design$fun,1,4)!="mrss"){
            idxx <- match(pars[idx[!is.na(idx)]], names(design$par))
            xseq <- seq(max(3, design$par[[idxx]]-5), design$par[[idxx]]+20, 1)
          }else{
            idxxx <- match(pars[idx[!is.na(idx)]],colnames(design$round.mrss))
            xseq <- seq(max(3, design$round.mrss[idxxx]-5), design$round.mrss[idxxx]+20, 1)
          }
        }else{
          xseq <- mrss.seq
        }
      }else if(x=="power"){
        if(is.null(power.seq)){
          if(!is.null(design$par$power)){
            xseq <- seq(max(.20, design$par$power-.50), min(.99, design$par$power+.20), .01)
          }else{
            xseq <- seq(max(.20, design$power-.50), min(.99, design$power+.20), .01)
          }
        }else{
          xseq <- power.seq
        }
      }
    }
    if(substr(design$fun,1,5)=="power" & "power"%in%pars){
      if(!is.null(power.seq)){
        stop("Sequence for 'power' is not allowed.")
      }
      y <- "power"
      x <- pars[pars!="power"][1]
      idx <- match(c("n","J","K","L"),pars)
      if(any(!is.na(idx))){
        if(is.null(mrss.seq)){
          if(substr(design$fun,1,4)!="mrss"){
            idxx <- match(pars[idx[!is.na(idx)]], names(design$par))
            xseq <- seq(max(3, design$par[[idxx]]-5), design$par[[idxx]]+20, 1)
          }else{
            idxxx <- match(pars[idx[!is.na(idx)]],colnames(design$round.mrss))
            xseq <- seq(max(3, design$round.mrss[idxxx]-5), design$round.mrss[idxxx]+20, 1)
          }
        }else{
          xseq <- mrss.seq
        }
      }else if(x=="mdes"){
        if(is.null(mdes.seq)){
          if(!is.null(design$par$mdes)){
            xseq <- seq(max(.05, design$par$mdes-.50), design$par$mdes+.20, .01)
          }else{
            xseq <- seq(max(.05, design$mdes[1]-.50), design$mdes[1]+.20, .01)
          }
        }else{
          xseq <- mdes.seq
        }
      }
    }

    # indexes
    idx <- match(x,names(design$par))
    x0 <- which(xseq==round(design$par[[idx]], digits=2))
    idy <- match(y,names(design))

    # plot
    yout <- matrix(NA,1,length(xseq))
    for(i in 1:length(xseq)){
      par <- design$par
      par[[idx]] <- xseq[i]
      yout[i] <- do.call(design$fun, par)[[idy]][1]
    }
    plot(xseq, yout,
         type="l",
         xlab=x, ylab=y,
         xlim=range(xseq),
         ylim=range(yout),
         ...)
    points(xseq[x0], yout[x0], col='black', pch=21, bg = "red", cex=1.5)
  }else{
    stop("Specify two or three parameters.")
  }

}

# examples
# design1 <- mdes.bcra4f3(rho3=.15, rho2=.15, n=10, J=4, L=27, K=4)
# design2 <- power.bcra4f3(rho3=.15, rho2=.15, n=10, J=4, L=27, K=4)
# plot.pars(design1)
# plot.pars(design2)

