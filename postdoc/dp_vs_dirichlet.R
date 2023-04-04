library(gtools)

nsamp <- 1E4
K <- 7
a <- 7

stick_breaking <- t(replicate(nsamp, {
  b <- rbeta(K, 1, a)
  sort(c(b[1], sapply(2:K, function(i) b[i] * prod(1 - b[1:(i-1)]))), decreasing = T)
}))
their_thing <- t(apply(gtools::rdirichlet(nsamp, rep(a/K, K)), 1, sort, decreasing = T))

par(mfrow = c(2,1))
hist(stick_breaking[,1])
hist(their_thing[,1])


mean(gtools::rdirichlet(nsamp, rep(a/K, K))[,1])

p <- 0:1000/1000
d1 <- extraDistr::dkumar(p, 1, 5)
d2 <- dbeta(p, 1, 5)

plot(d1, d2)
max(abs(d1 - d2))
