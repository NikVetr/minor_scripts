library(rethinking)

#5M4
data(WaffleDivorce)
WaffleDivorce <- WaffleDivorce[order(as.character(WaffleDivorce$Location)),]

LDS <- read.csv("LDS.csv")
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

m5M4 <- quap( 
  alist(
  D ~ dnorm(mu, sigma),
  mu <- a + M*bM + A*bA + LDS*bLDS,
  c(bM, bA, bLDS) ~ dnorm(0, 0.5),
  a ~ dnorm(0,0.3),
  sigma ~ dexp(1)
), data = d)
precis(m5M4)
plot(coeftab(m5M4))

#5H1
data(foxes)
str(foxes)

d <- list()
d$g <- foxes$group
d$bw <- standardize(foxes$weight)
d$a <- standardize(foxes$area)
d$gs <- standardize(foxes$groupsize)

m5H1A <- quap(
  alist(
    bw ~ dnorm(mu, sig),
    mu <- i + a*bA,
    a ~ dnorm(0,0.3),
    c(i, bA) ~ dnorm(0,0.3),
    sig ~ dexp(1)
  ),
  data = d
)

m5H1B <- quap(
  alist(
    bw ~ dnorm(mu, sig),
    mu <- i + gs*bGS,
    c(i, bGS) ~ dnorm(0,0.3),
    sig ~ dexp(1)
  ),
  data = d
)

plot(coeftab(m5H1A, m5H1B))

#mA
sim_a <- 1:100 / 25  - 2
muA <- link(m5H1A, data = list(a = sim_a))
muA_mu <- apply(X = muA, MARGIN = 2, mean)
muA_PI <- apply(X = muA, MARGIN = 2, PI)
plot(sim_a, muA_mu, type = "l", lwd = 2, col = 2, ylim = c(-2,2))
shade(muA_PI, sim_a)

#mB
sim_gs <- 1:100 / 25  - 2
muB <- link(m5H1B, data = list(gs = sim_gs))
muB_mu <- apply(X = muB, MARGIN = 2, mean)
muB_PI <- apply(X = muB, MARGIN = 2, PI)
plot(sim_gs, muB_mu, type = "l", lwd = 2, col = 2, ylim = c(-2,2))
shade(muB_PI, sim_gs)

#5H2
m5H2 <- quap(
  alist(
    bw ~ dnorm(mu, sig),
    mu <- i + bA*a + bGS*gs,
    c(bA, bGS) ~ dnorm(0, 0.5),
    i ~ dnorm(0,0.3),
    sig ~ dexp(1)
  ),
  data = d
)
plot(coeftab(m5H2))
precis(m5H2)
sim_gs <- sim(m5H2, data = list(a = rep(0, 100), gs = 1:100 / 25  - 2))
sim_a <- sim(m5H2, data = list(gs = rep(0, 100), a = 1:100 / 25  - 2))

sim_gs_mu <- apply(X = sim_gs, MARGIN = 2, mean)
sim_gs_PI <- apply(X = sim_gs, MARGIN = 2, PI)
plot(1:100 / 25  - 2, sim_gs_mu, type = "l", lwd = 2, col = 2, ylim = c(-2,2))
shade(sim_gs_PI, 1:100 / 25  - 2)


sim_a_mu <- apply(X = sim_a, MARGIN = 2, mean)
sim_a_PI <- apply(X = sim_a, MARGIN = 2, PI)
plot(1:100 / 25  - 2, sim_a_mu, type = "l", lwd = 2, col = 2, ylim = c(-2,2))
shade(sim_a_PI, 1:100 / 25  - 2)

plot(d$a, d$gs)

#5H3

d$f <- standardize(foxes$avgfood)

m5H3A <- quap(
  alist(
    bw ~ dnorm(mu, sig),
    mu <- i + gs*bGS + f*bF,
    c(i, bGS, bF) ~ dnorm(0,0.3),
    sig ~ dexp(1)
  ),
  data = d
)

m5H3B <- quap(
  alist(
    bw ~ dnorm(mu, sig),
    mu <- i + gs*bGS + f*bF + bA*a,
    c(i, bGS, bF, bA) ~ dnorm(0,0.3),
    sig ~ dexp(1)
  ),
  data = d
)

plot(coeftab(m5H1A, m5H1B, m5H2, m5H3A, m5H3B))
plot(coeftab(m5H3A, m5H3B))
pairs(d)
pairs(m5H3B)
#neither, both associated with outcome and eachother, together, parameters are nonidentifiable
