n <- 5
p <- 2
r <- 0.5
R <- diag(n) * (1-r) + r
L <- t(chol(R))
detR <- det(R)

ttsim <- function(n, p, L, detR){
  x <- L %*% matrix(rnorm(n * p), n, p)
  mu <- apply(x,2,mean)
  dmu <- diff(mu)
  tval <- dmu / sqrt(sum(apply(x, 2, function(i) var(i) / length(i))))
  c(myp = 1 - abs(0.5 - pt(tval, df = 0 + detR * (n * p - 4))) * 2,
    ttp = t.test(x[,1], x[,2])$p.value)
}

nrep <- 1E4
out <- data.frame(t(replicate(nrep, ttsim(n, p, L, detR))))

par(mfrow = c(3,1))
hist(out$myp, xlim = c(0,1), breaks = 0:10/10)
hist(out$ttp, xlim = c(0,1), breaks = 0:10/10)
plot(out[order(out$myp),], type = "l")
cor(out)
