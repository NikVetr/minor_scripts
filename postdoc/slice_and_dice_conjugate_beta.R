x <- matrix(rbinom(4, 100, rbeta(4, 3, 3)), 2, 2)
n <- 1E6
pp1 <- rbeta(n, 1 + x[1,1], 1 + x[1,2]) - rbeta(n, 1 + x[2,1], 1 + x[2,2])
pp2 <- rbeta(n, 1 + x[1,1], 1 + x[2,1]) - rbeta(n, 1 + x[1,2], 1 + x[2,2])
both <- c(pp1, pp2)
breaks <- seq(from = min(both) - diff(range(both)) / 100, to =  max(both) + diff(range(both)) / 100, length.out = 20)

hist(pp1, col = adjustcolor(2,0.5), breaks = breaks)
hist(pp2, col = adjustcolor(3,0.5), breaks = breaks, add = T)


mean(pp1)
mean(pp2)
mean(pp1 > 0)
mean(pp2 > 0)

qs <- 0:100 / 100
plot(quantile(pp1, qs), quantile(pp2, qs), type = "l", lwd = 2)
abline(0,1,col = 3, lty = 2)
