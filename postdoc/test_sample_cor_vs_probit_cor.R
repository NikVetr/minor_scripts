rij <- -0.8
R <- matrix(c(1,rij,rij,1), 2, 2)
n <- 1E3
x <- matrix(rnorm(n*2), 2, n)
x <- t(t(chol(R)) %*% x)
plot(x, pch = 19, col = adjustcolor(1,0.5))

thresh <- 2
y <- (x > thresh) + 0

cor(x)
cor(y)

table(apply(y, 1, paste0, collapse = "-"))

z <- y[,1]
zbar <- mean(z)
sum(log(1 - apply(cbind(zbar, abs(z - zbar)), 1, max)))

z_est <- 0:1000/1000
ll <- sapply(z_est, function(zi) sum(log(1 - apply(cbind(zi, abs(z - zi)), 1, max))))
plot(z_est, ll, type = "l")
z_est[which.max(ll)]
zbar
