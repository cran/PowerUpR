## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(PowerUpR)

## ----  message=FALSE, eval=FALSE-----------------------------------------
#  install.packages("PowerUpR")

## ---- message=FALSE------------------------------------------------------
mdes.cra3r3(power=.80, rho2=.06, rho3=.18,
            g3=1, r21=.55, r22=.50, r23=.45,
            p=.40, n=10, J=2, K=83)

## ---- message=FALSE------------------------------------------------------
power.cra3r3(es=.23, rho2=.06, rho3=.18,
             g3=1, r21=.55, r22=.50, r23=.45,
             p=.40, n=10, J=2, K=83)


