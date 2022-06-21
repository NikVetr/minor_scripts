set.seed(5)

n <- 100
x <- rnorm(n)
alpha <- 2
beta <- 0
epsilon_sd <- 1

y <- alpha + x * beta + rnorm(n, sd = epsilon_sd)
fit1 <- lm(y~x)$coefficients
fit2 <- lm(x~y)$coefficients
fit2 <- c(-fit2[1] / fit2[2], 1 / fit2[2])
v <- prcomp(cbind(x,y))$rotation
fit3 <- v[2,1]/v[1,1]
fit3 <- c(mean(y) - fit3 * mean(x), fit3)

plot(x, y)
abline(fit1[1], fit1[2], col = 2, lwd = 2)
abline(fit2[1], fit2[2], col = 4, lwd = 2)
abline(fit3[1], fit3[2], col = "purple", lwd = 2)

fits <- sapply(list(fit1, fit2, fit3), function(i) 
  paste0("y = ", round(i[1], 2), " + ", round(i[2], 2), " * x"))

legend(x = "topleft", legend = c("y ~ x (ols)", "x ~ y (ols, inverse)", "y ~ x (tls)"),
       lwd = 2, col = c(2,4,"purple"))
legend(x = "bottomright", legend = fits,
       lwd = 2, col = c(2,4,"purple"))

