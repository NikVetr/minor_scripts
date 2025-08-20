library(pROC)
set.seed(1)

#simulate training data
n <- 1E3
x_sd <- 1
e_sd <- 0.8
b <- 1
x <- rnorm(n, 0, x_sd)
e <- rnorm(n, 0, e_sd)
y <- x * b + e

#fit model
fit <- lm(y ~ 1 + x)
sfit <- summary(fit)
sfit

#simulate new data
n2 <- 1E5
x2 <- rnorm(n2, 0, x_sd)
e2 <- rnorm(n2, 0, e_sd)
y2 <- x2 * b + e2

#make predictions on new data
mu_hat <- predict(fit, newdata = data.frame(x = x2))  
sigma_hat <- summary(fit)$sigma

#compare to threshold
qthresh <- 0.01
thresh <- quantile(y, qthresh)
p_above_thresh <- 1 - pnorm(thresh, mean = mu_hat, sd = sigma_hat)
roc_obj <- roc(y2 > thresh, p_above_thresh, quiet = TRUE)
df <- data.frame(
  threshold = roc_obj$thresholds,
  fpr       = 1 - roc_obj$specificities,
  fnr       = 1 - roc_obj$sensitivities
)
target_fpr <- 0.05
idx <- which.min(abs(df$fpr - target_fpr))

#plot result
plot(df$fpr, df$fnr, type = "l",
     xlab = "False-positive rate", ylab = "False-negative rate",
     xlim = c(0,1), ylim = c(0,1), asp = 1)

#label the point corresponding to an fpr of 0.05
points(df$fpr[idx], df$fnr[idx], pch = 19, col = 2)
text(df$fpr[idx], df$fnr[idx],
     labels = sprintf("threshold = %.3f", df$threshold[idx]),
     pos = 4, col = 2)
segments(x0 = target_fpr, x1 = target_fpr, y0 = par("usr")[3], y1 = df$fnr[idx], col = 2, lty = 2)
segments(x0 = par("usr")[1], x1 = target_fpr, y0 = df$fnr[idx], y1 = df$fnr[idx], col = 2, lty = 2)

pred_lab <- as.integer(p_above_thresh > df$threshold[idx])
true_lab <- as.integer(y2 > thresh)
cm <- table(Predicted = pred_lab,
            Actual    = true_lab)
epitools::oddsratio(cm)$measure

