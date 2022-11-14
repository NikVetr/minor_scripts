rn <- 1000

out <- data.frame(t(replicate(rn, {
  n <- 10
  sd_x <- 2
  x <- rnorm(n, sd = sd_x)
  sd_dx <- 0.5
  mean_dx <- 0.25
  dx <- rnorm(n, mean = mean_dx, sd = sd_dx)
  y <- x + dx
  
  #add noise to out observed x and y
  x <- x + rnorm(n, sd = 1)
  y <- y + rnorm(n, sd = 1)
  
  #run both tests
  c(two.sample = t.test(x=x, y=y)$p.value, one.sample = t.test(y-x)$p.value)
})))

alpha <- 0.05
cols <- rep(1, rn)
cols[out$one.sample < alpha] <- 2
cols[out$two.sample < alpha] <- 4
cols[out$one.sample < alpha & out$two.sample < alpha] <- 6
plot(out, xlim=c(0,1), ylim=c(0,1), main = "p-value comparison for 1 vs 2 sample tests",
     pch = 19, col = adjustcolor(cols, 0.5)); 
abline(h = alpha, col = 2, lty = 2); abline(v = alpha, col = 4, lty = 2)

hist(out$two.sample)
hist(out$one.sample)
