#load packages
library(rethinking)
library(mcmcplots)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

d <- read.csv("C:\\Users\\Nikolai\\Dropbox\\R Files\\Tomasik_Data.csv")
d <- d[which(d$RequesterFeedback == ""), 23:33]
str(d)
d$presRain <- coerce_index(d$Answer.preserve.rainforest) - 1
d$animSuff <- coerce_index(d$Answer.animal.suffering) - 1
d$male <- coerce_index(d$Answer.Gender) - 1
d$ageC <- d$Answer.Age - mean(d$Answer.Age) #center to ease later interpretation
d$labUniv <- coerce_index(d$Answer.lab.universes) - 1
d$space <- coerce_index(d$Answer.spread.life.to.other.planets) - 1

#Preserve the Rainforest

m0 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI, 
    pI ~ dnorm(0,1)
  ) ,
  data= list(presRain = d$presRain)
)
precis(m0)

m1 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS, 
    c(pI, pAS) ~ dnorm(0,1)
  ) ,
  data= list(presRain = d$presRain, animSuff = d$animSuff)
)
precis(m1)

m2 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM, 
    c(pI, pAS, pM) ~ dnorm(0,1)
  ) ,
  data= list(presRain = d$presRain, animSuff = d$animSuff, male <- d$male)
)
precis(m2)

m3 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + animSuff * male * pAMS, 
    c(pI, pAS, pM, pAMS) ~ dnorm(0,1)
  ) ,
  data= list(presRain = d$presRain, animSuff = d$animSuff, male = d$male)
)
precis(m3)

m4 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + animSuff * male * pAMS + ageC * pA, 
    c(pI, pAS, pM, pAMS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(presRain = d$presRain, animSuff = d$animSuff, male = d$male, ageC = d$ageC)
)
precis(m4)

m5 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + ageC * pA, 
    c(pI, pAS, pM) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(presRain = d$presRain, animSuff = d$animSuff, male = d$male, ageC = d$ageC)
)
precis(m5)

m6 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + male*pM, 
    c(pI, pM) ~ dnorm(0,1)
  ) ,
  data= list(presRain = d$presRain, male = d$male)
)
precis(m6)

m7 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA, 
    pI ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(presRain = d$presRain, ageC = d$ageC)
)
precis(m7)

m8 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + animSuff*pAS, 
    c(pI, pAS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(presRain = d$presRain, ageC = d$ageC, animSuff = d$animSuff)
)
precis(m8)

m9 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + animSuff*pAS + pAAS*animSuff*ageC, 
    c(pI, pAS, pAAS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(presRain = d$presRain, ageC = d$ageC, animSuff = d$animSuff)
)
precis(m9)


m10 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + animSuff*pAS + pAAS*animSuff*ageC + male*pM, 
    c(pI, pAS, pAAS, pM) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(presRain = d$presRain, ageC = d$ageC, animSuff = d$animSuff, male = d$male)
)
precis(m10)

m11 <- map( 
  alist(
    presRain ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + animSuff*pAS + pAAS*animSuff*ageC + male*pM + pMAS*male*animSuff , 
    c(pI, pAS, pAAS, pM, pMAS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(presRain = d$presRain, ageC = d$ageC, animSuff = d$animSuff, male = d$male)
)
precis(m11)


compare(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11)

##Create Lab Universes##

m0 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI, 
    pI ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv)
)
precis(m0)

m1 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS, 
    c(pI, pAS) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff)
)
precis(m1)

m2 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM, 
    c(pI, pAS, pM) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, male <- d$male)
)
precis(m2)

m3 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + animSuff * male * pAMS, 
    c(pI, pAS, pM, pAMS) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, male = d$male)
)
precis(m3)

m4 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + animSuff * male * pAMS + ageC * pA, 
    c(pI, pAS, pM, pAMS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, male = d$male, ageC = d$ageC)
)
precis(m4)

m5 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + ageC * pA, 
    c(pI, pAS, pM) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, male = d$male, ageC = d$ageC)
)
precis(m5)

m6 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + male*pM, 
    c(pI, pM) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, male = d$male)
)
precis(m6)

m7 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA, 
    pI ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, ageC = d$ageC)
)
precis(m7)

m8 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + pS*Space, 
    c(pI, pS) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, Space = d$space)
)
precis(m8)

m9 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS  + pS*Space, 
    c(pI, pAS, pS) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, Space = d$space)
)
precis(m9)

m10 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM +  + pS*Space, 
    c(pI, pAS, pM, pS) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, male <- d$male, Space = d$space)
)
precis(m10)

m11 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + animSuff * male * pAMS +  + pS*Space, 
    c(pI, pAS, pM, pAMS, pS) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, male = d$male, Space = d$space)
)
precis(m11)

m12 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + animSuff * male * pAMS + ageC * pA +  + pS*Space, 
    c(pI, pAS, pM, pAMS, pS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, male = d$male, ageC = d$ageC, Space = d$space)
)
precis(m12)

m13 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + animSuff*pAS + male*pM + ageC * pA +  + pS*Space, 
    c(pI, pAS, pM, pS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, animSuff = d$animSuff, male = d$male, ageC = d$ageC, Space = d$space)
)
precis(m13)

m14 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + male*pM +  + pS*Space, 
    c(pI, pM, pS) ~ dnorm(0,1)
  ) ,
  data= list(labUniv = d$labUniv, male = d$male, Space = d$space)
)
precis(m14)

m15 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + pS*Space, 
    c(pI, pS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, ageC = d$ageC, Space = d$space)
)
precis(m15)

m16 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + animSuff*pAS, 
    c(pI, pAS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, ageC = d$ageC, animSuff = d$animSuff)
)
precis(m16)

m17 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + animSuff*pAS + pAAS*animSuff*ageC, 
    c(pI, pAS, pAAS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, ageC = d$ageC, animSuff = d$animSuff)
)
precis(m17)


m18 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + animSuff*pAS + pAAS*animSuff*ageC + male*pM, 
    c(pI, pAS, pAAS, pM) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, ageC = d$ageC, animSuff = d$animSuff, male = d$male)
)
precis(m18)

m19 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + ageC * pA + animSuff*pAS + pAAS*animSuff*ageC + male*pM + pMAS*male*animSuff , 
    c(pI, pAS, pAAS, pM, pMAS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, ageC = d$ageC, animSuff = d$animSuff, male = d$male)
)
precis(m19)

compare(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19)

m20 <- map( 
  alist(
    labUniv ~ dbinom(1, p),
    logit(p) <- pI + male*pM + ageC * pA +  + pS*Space, 
    c(pI, pM, pS) ~ dnorm(0,1),
    pA ~ dnorm(0, 0.5)
  ) ,
  data= list(labUniv = d$labUniv, male = d$male, ageC = d$ageC, Space = d$space)
)
precis(m20)

compare(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20)

