n <- 100
x <- rnorm(n)
a <- 2
b <- 1.5
y <- a + b * x
e <- lm(rnorm(n) ~ y)$residuals
e <- e / sd(e)
e <- rnorm(n)
ye <- y + e

fit <- lm(ye~1+x)
rye <- ye - cov(x,ye)/var(x) * x

summary(fit)$coefficients
mean(rye)
vcov(fit)
sd(rye) / sqrt(n)
