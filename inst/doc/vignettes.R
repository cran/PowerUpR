## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- message=FALSE, warning=FALSE---------------------------------------
library(PowerUpR)

## ---- message=FALSE------------------------------------------------------
design1 <- mdes.cra3r3(power=.80, rho2=.06, rho3=.18,
                       g3=1, R12=.55, R22=.50, R32=.45,
                       P=.40, n=10, J=2, K=83)
print(design1$mdes)

## ---- message=FALSE------------------------------------------------------
design2 <- power.cra3r3(mdes=.23, rho2=.06, rho3=.18,
                        g3=1, R12=.55, R22=.50, R32=.45,
                        P=.40, n=10, J=2, K=83)
print(design2$power)

## ---- message=FALSE------------------------------------------------------
design3 <- mrss.cra3r3(power=.80, mdes=.23, rho2=.06, rho3=.18,
                       g3=1, R12=.55, R22=.50, R32=.45,
                       P=.40, n=10, J=2)
print(design3$round.mrss)

## ---- message=FALSE------------------------------------------------------
design4 <- mrss.cra3r3(power=.80, mdes=.23, rho2=.06, rho3=.18,
                       g3=1, R12=.55, R22=.50, R32=.45,
                       P=.40, n=10, K=83)
print(design4$round.mrss)

## ---- message=FALSE------------------------------------------------------
design5 <- mrss.cra3r3(power=.80, mdes=.23, rho2=.06, rho3=.18,
                       g3=1, R12=.55, R22=.50, R32=.45,
                       P=.40, J=2, K=83)
print(design5$round.mrss)

## ---- message=FALSE------------------------------------------------------
design6 <- optimal.cra3r3(cn=10, cJ=200, cK=500, constrain="power", power=.80, gm=5, ncase=20,
                          mdes=.23, rho2=.06, rho3=.18,
                          g3=1, R12=.55, R22=.50, R32=.45,
                          P=.40)
print(design6$round.optim)
print(design6$integer.optim)

## ---- message=FALSE------------------------------------------------------
design7 <- optimal.cra3r3(cn=10, cJ=200, cK=500, constrain="power", power=.80, J=2, 
                          mdes=.23, rho2=.06, rho3=.18,
                          g3=1, R12=.55, R22=.50, R32=.45,
                          P=.40)
print(design7$round.optim)
print(design7$integer.optim)

## ---- message=FALSE------------------------------------------------------
design8 <- optimal.cra3r3(cn=10, cJ=200, cK=500, constrain="mdes", mdes=.20, J=2,
                          power=.80, rho2=.06, rho3=.18,
                          g3=1, R12=.55, R22=.50, R32=.45,
                          P=.40)
print(design8$round.optim)
print(design8$integer.optim)

## ---- message=FALSE------------------------------------------------------
design9 <- optimal.cra3r3(cn=10, cJ=200, cK=500, constrain="cost", cost=130000, J=2,
                          power=.80, mdes=.20, rho2=.06, rho3=.18,
                          g3=1, R12=.55, R22=.50, R32=.45,
                          P=.40)
print(design9$round.optim)
print(design9$integer.optim)

## ---- message=FALSE------------------------------------------------------
design10 <- optimal.cra3r3(cn=10, cJ=200, cK=500, constrain="cost", cost=130000, n=20, J=2,
                           power=.80, mdes=.20, rho2=.06, rho3=.18,
                           g3=1, R12=.55, R22=.50, R32=.45,
                           P=.40)
print(design10$round.optim)
print(design10$integer.optim)

## ---- message=FALSE------------------------------------------------------
## MDES - MRSS(K)
## Red point indicates the design location
plot(design10, pars=c("mdes","K"))

## ---- message=FALSE------------------------------------------------------
## Power - MRSS(K)
## Red point indicates the design location
plot(design10, pars=c("power","K"))

## ---- message=FALSE------------------------------------------------------
## Perspective plot
## Red point indicates the design location
plot(design10, type="p")

## ---- message=FALSE------------------------------------------------------
## Contour plot
## Red point indicates the design location
plot(design10, type="c")

## ---- message=FALSE------------------------------------------------------
## confidence intervals for MDES
design11 <- optimal.to.mdes(design10)
print(design11$mdes)
## percentile equivalance
mdes.to.pctl(design11)


