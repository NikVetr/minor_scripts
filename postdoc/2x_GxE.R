logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x) / (1 + exp(x))

n1 <- 5000
n2 <- 5000
p <- 1000
f1 <- rbeta(p, 5, 20)
f2 <- invlogit(logit(f1) - 0.5)

x1 <- sapply(f1, function(prob) rbinom(n1, 2, prob))
x2 <- sapply(f2, function(prob) rbinom(n2, 2, prob))
ind <- c(rep(1, n1), rep(2, n2))

x <- rbind(x1, x2)

h2_tot <- 0.3
b <- t(t(rnorm(p)))

yg <- c(x %*% b)
y <- yg + rnorm(length(yg), sd = sqrt(var(yg) * (1 - h2_tot) / h2_tot))

y1 <- y[ind == 1]
y2 <- y[ind == 2]

fits <- sapply(1:p, function(pi) {
  fit1 <- lm(y1 ~ x1[,pi])
  fit2 <- lm(y2 ~ x2[,pi])
  
  c(x1_coef = summary(fit1)$coefficients[2,c(1,4)],
    x2_coef = summary(fit2)$coefficients[2,c(1,4)])
})

fits <- data.frame(t(fits))

alpha <- 0.1
pop1_sig <- which(p.adjust(fits$x1_coef.Pr...t.., method = "BH") < alpha)
pop2_sig <- which(p.adjust(fits$x2_coef.Pr...t.., method = "BH") < alpha)

both_sig <- intersect(pop1_sig, pop2_sig)
# both_sig <- which(fits$x1_coef.Pr...t.. < alpha & fits$x2_coef.Pr...t.. < alpha)

plot(fits$x1_coef.Estimate[both_sig], fits$x2_coef.Estimate[both_sig], col = adjustcolor(1, 0.5), pch = 19)
abline(0,1, col = "grey50", lty = 2, lwd = 2)
linefit <- summary(lm(fits$x2_coef.Estimate[both_sig] ~ fits$x1_coef.Estimate[both_sig]))$coefficients
(1 - linefit[2,1]) / linefit[2,2]
abline(a = linefit[1,1], b = linefit[2,1], col = "red", lty = 3, lwd = 2)

v <- prcomp(cbind(fits$x1_coef.Estimate[both_sig],fits$x2_coef.Estimate[both_sig]))$rotation
tls_fit <- v[2,1]/v[1,1]
tls_fit <- c(mean(fits$x2_coef.Estimate[both_sig]) - tls_fit * mean(fits$x1_coef.Estimate[both_sig]), tls_fit)
abline(tls_fit[1], tls_fit[2], col = "purple", lwd = 2, lty = 3)

tls_boots <- t(replicate(1000, {
    boot_inds <- sample(1:length(fits$x1_coef.Estimate[both_sig]), length(fits$x1_coef.Estimate[both_sig]), T)
    v <- prcomp(cbind(fits$x1_coef.Estimate[both_sig][boot_inds],fits$x2_coef.Estimate[both_sig][boot_inds]))$rotation
    tls_fit_boot <- v[2,1]/v[1,1]
    tls_fit_boot <- c(mean(fits$x2_coef.Estimate[both_sig][boot_inds]) - tls_fit_boot * mean(fits$x1_coef.Estimate[both_sig][boot_inds]), tls_fit_boot)
    tls_fit_boot
}))

hist(tls_boots[,2], col = "purple")
