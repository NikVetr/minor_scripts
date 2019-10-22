g1 <- rnorm(50, mean = 80, sd = 15); g1[g1>100] <- 100; g1[g1<0] <- 0
g2 <- rnorm(50, mean = 20, sd = 15); g2[g2>100] <- 100; g2[g2<0] <- 0
b2 <- rnorm(50, mean = 80, sd = 15); b1[b1>100] <- 100; b2[b2<0] <- 0
b1 <- rnorm(50, mean = 20, sd = 15); b2[b2>100] <- 100; b1[b1<0] <- 0

plot(c(b1,g1), c(b2,g2), xlab = "Hockey-Playing Ability", ylab = "French-Speaking Ability")
abline(lm(c(b2,g2) ~ c(b1,g1)), lwd = 2)
fit <- (summary(lm(c(b2,g2) ~ c(b1,g1))))
fit$r.squared
p.val <- signif(fit$coefficients[2,4]* 1e100, 3) / 1e100
text(x = 80, y = 95, paste0("p-value = ", p.val), cex = 3)

plot(c(b1,g1), c(b2,g2), xlab = "Hockey-Playing Ability", ylab = "French-Speaking Ability")
abline(lm(c(b2,g2) ~ c(b1,g1)), lwd = 2)
fit <- (summary(lm(c(b2,g2) ~ c(b1,g1))))
fit$r.squared
p.val <- signif(fit$coefficients[2,4]* 1e100, 3) / 1e100
text(x = 80, y = 95, paste0("p-value = ", p.val), cex = 3)
intercept <- lm(c(b2,g2) ~ c(b1,g1))$coefficients[[1]]
slope <- lm(c(b2,g2) ~ c(b1,g1))$coefficients[[2]]
x <- c(b1,g1)
y <- c(b2,g2)
for(i in 1:length(c(b1,g1))){
  y.pred <- intercept + x[i]*slope
  segments(x0 = x[i], y0 = y[i], x1 = x[i], y1 = y.pred, lty = 2)
}

plot(-40:40/10, dnorm(-40:40/10), type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = 3, main = "Normal Distribution", cex.main = 3)

plot(c(b1,g1), c(b2,g2), xlab = "Hockey-Playing Ability", ylab = "French-Speaking Ability", col = c(rep("red", 50), rep("blue", 50)), pch = 16)
abline(lm(c(b2,g2) ~ c(b1,g1)), lwd = 2)
fit <- (summary(lm(c(b2,g2) ~ c(b1,g1))))
p.val <- signif(fit$coefficients[2,4]* 1e100, 3) / 1e100
text(x = 80, y = 95, paste0("p-value = ", p.val), cex = 3)

plot(c(b1,g1), c(b2,g2), xlab = "Hockey-Playing Ability", ylab = "French-Speaking Ability", col = c(rep("red", 50), rep("blue", 50)), pch = 16)
abline(lm(c(b2,g2) ~ c(b1,g1)), lwd = 1, lty = 2)
intercept.1 <- lm(c(b2) ~ c(b1))$coefficients[[1]]
slope.1 <- lm(c(b2) ~ c(b1))$coefficients[[2]]
intercept.2 <- lm(c(g2) ~ c(g1))$coefficients[[1]]
slope.2 <- lm(c(g2) ~ c(g1))$coefficients[[2]]
segments(x0 = 0, x1 = 50, y0 = intercept.1, y1 = intercept.1 + 50*slope.1, lwd = 2)
segments(x0 = 50, x1 = 100, y0 = intercept.2, y1 = intercept.2 + 50*slope.2, lwd = 2)
#this isn't right but gets the point across
fit <- summary(lm(c(b1,g1) ~ c(b2-mean(b2),g2-mean(g2))))
p.val <- signif(fit$coefficients[2,4]* 1e100, 3) / 1e100
text(x = 80, y = 95, paste0("p-value = ", p.val), cex = 3)



