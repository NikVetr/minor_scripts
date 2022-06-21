library(rethinking)

data(WaffleDivorce)
?WaffleDivorce
WaffleDivorce <- WaffleDivorce[order(as.character(WaffleDivorce$Location)),]

library(dagitty)
dagW <- dagitty( "dag {
A -> D
A -> M -> D
A <- S -> M
S -> W -> D
}")

coordinates(dagW) <- list(x = c(A=rnorm(1), D=rnorm(1), M=rnorm(1), S=rnorm(1), W =rnorm(1)) ,
                         y = c(A=rnorm(1), D=rnorm(1), M=rnorm(1), S=rnorm(1), W =rnorm(1)))
plot(dagW)


LDS <- read.csv("data/LDS.csv")
LDS <- LDS[order(as.character(LDS$State)),]; LDS <- LDS[-which(LDS$State == "Nevada"),]
LDS$Membership <- as.character(LDS$Membership)
LDS$Membership <- as.numeric(gsub(pattern = ",", replacement = "", x = LDS$Membership))
LDS$Population <- as.character(LDS$Population)
LDS$Population <- as.numeric(gsub(pattern = ",", replacement = "", x = LDS$Population))


d <- list()
d$A <- standardize( WaffleDivorce$MedianAgeMarriage )
d$D <- standardize( WaffleDivorce$Divorce )
d$M <- standardize( WaffleDivorce$Marriage )
d$LDS <- standardize(LDS$Membership / LDS$Population)
d$S <- WaffleDivorce$South
d$W <- standardize(WaffleDivorce$WaffleHouses / WaffleDivorce$Population)

panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 4 * abs(cor(x, y)))
}

pairs(d, upper.panel = panel.cor)

m6.1 <- quap( 
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + W*bW + S*bS,
    c(bW, bS) ~ dnorm(0, 0.5),
    a ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.1)
plot(coeftab(m6.1))

d$PS <- standardize(WaffleDivorce$PropSlaves1860)
m6.2 <- quap( 
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + W*bW + S*bS + PS*bPS,
    c(bW, bS, bPS) ~ dnorm(0, 0.5),
    a ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.2)
plot(coeftab(m6.2))


m6.3 <- quap( 
  alist(
    A ~ dnorm(mu, sigma),
    mu <- a + W*bW + S*bS,
    c(bW, bS) ~ dnorm(0, 0.5),
    a ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.3)
plot(coeftab(m6.3))
#W and A conditionally independent

m6.4 <- quap( 
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + W*bW + S*bS,
    c(bW, bS) ~ dnorm(0, 0.5),
    a ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.4)
plot(coeftab(m6.4))

m6.5 <- quap( 
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + W*bW + S*bS + PS*bPS + LDS*bLDS,
    c(bW, bS, bPS, bLDS) ~ dnorm(0, 0.5),
    a ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data = d)
precis(m6.5)
plot(coeftab(m6.5))
