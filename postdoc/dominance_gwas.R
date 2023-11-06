n <- 1E3
x <- rbinom(n, 2, 0.75)
one <- (x > 0) * 1
two <- (x > 1) * 1
y <- one + two * 0.5 + rnorm(n) * 2

coefficients(summary(lm(y ~ one + two)))
coefficients(summary(lm(y ~ two)))
coefficients(summary(lm(y ~ x)))
