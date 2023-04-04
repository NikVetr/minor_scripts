n <- 100
mus <- c(-10, 10)
sds <- c(1,1)
d <- data.frame(y = c(rnorm(n, mus[1], sds[1]), rnorm(n, mus[2], sds[2])), x = as.factor(rep(1:2, each = n)))
d$y[1] <- 1000
d$y <- lm(d$y ~ d$x)$resid

plot(d$y, d$x, pch = 19, col = d$x, xlim = range(d$y[-1]))

mean(resids[d$x == 1])
mean(resids[d$x == 2])

table(d$x[resids > 0])
table(d$x[resids < 0])

library(e1071)
tune.res <- tune(svm, x ~ y, data = d,  kernel = "linear", 
                 ranges = list(cost=10^(-3:3)))
best.cost <- as.numeric(tune.res$best.parameters)
fit <- svm(x ~ y, data = d, kernel = "linear", cost = best.cost)
preds <- predict(fit)
table(preds)

results <- caret::confusionMatrix(reference = d$x, 
                                  data = preds)
results
