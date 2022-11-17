# set.seed(1)
n <- 1E2 #number of individuals (eg humans) in the sample
m <- 5E2 #number of analytes / outcomes
p <- 50 #number of variables that systematically structure outcomes
k <- 15 #number of those variables that are known (and of interest)
n_PCs_to_use <- 29 #number of PCs to use
b <- matrix(rnorm(m*p), p, m) #true coefficients, effect of variables on each outcome
b[1:k,] <- 0 #the "null" hypothesis that those k variables are 0
# b[1:ceiling(k/2),] <- 0 #the "null" hypothesis that those k variables are 0
b[k+1,] <- 10 #make first unobserved covariate super strung
x <- matrix(rnorm(n*p), n, p) #values for each variable for each individual
y <- x %*% b #y observed without iid noise
e <- matrix(rnorm(n*m, sd = sqrt(p/2)), n, m) #residual error for y, with ~1/3rd variance indep
ye <- y + e #the y we observe, with residual noise

#incorporating no latent factors
fits_1 <- lapply(1:m, function(i){
  fit <- lm(ye[,i] ~ x[,1:k])
  return(list(
    coefs = summary(fit)$coefficients[-1,1],
    pvals = summary(fit)$coefficients[-1,4],
    resids = fit$resid
  ))
})
coefs_1 <- sapply(fits_1, function(x) x$coefs)
pvals_1 <- sapply(fits_1, function(x) x$pvals)
resids_1 <- sapply(fits_1, function(x) x$resids)

#then we do PCA on the residuals from the above
pca <- prcomp(resids_1)
# n_PCs_to_use <- sum(cumsum(prcomp(x)$sdev^2) / sum(prcomp(x)$sdev^2) < 0.9)
PCs_to_use <- 1:n_PCs_to_use
latent_vals <- pca$x[,PCs_to_use]

#show PCs are orthogonal to original variation
cor(cbind(x[,1:k], pca$x))[1:k, (k+1):p]
sum(cumsum(eigen(cor(resids_1))$values) / m < 0.9999) + 1
n - k - 1 #with an error term
p - k #with no error term

pc_corr_with_unobs_vars <- cor(cbind(latent_vals, x[,(k+1):p]))[1:n_PCs_to_use, (n_PCs_to_use+1):(n_PCs_to_use + p - k)]
spooky_latent_correlations <- t(apply(pc_corr_with_unobs_vars, 1, function(x){
  top_corr_ind <- which.max(abs(x))
  c(latent_var = top_corr_ind, 
    corr = max(abs(x)) * sign(x[top_corr_ind]),
    second_top_corr = x[order(abs(x), decreasing = T)[2]])
}))

#indeed these "latent factors" do covary with our unobserved variables
spooky_latent_correlations

#so let's include them in the regression
fits_2 <- lapply(1:m, function(i){
  fit <- lm(ye[,i] ~ x[,1:k] + latent_vals)
  return(list(
    coefs = summary(fit)$coefficients[2:(k+1),1],
    pvals = summary(fit)$coefficients[2:(k+1),4],
    resids = fit$resid
  ))
})
coefs_2 <- sapply(fits_2, function(x) x$coefs)
pvals_2 <- sapply(fits_2, function(x) x$pvals)
resids_2 <- sapply(fits_2, function(x) x$resids)

#compare to incl random noise
fits_3 <- lapply(1:m, function(i){
  fit <- lm(ye[,i] ~ x[,1:k] + matrix(rnorm(n*n_PCs_to_use), n, n_PCs_to_use))
  return(list(
    coefs = summary(fit)$coefficients[2:(k+1),1],
    pvals = summary(fit)$coefficients[2:(k+1),4],
    resids = fit$resid
  ))
})
coefs_3 <- sapply(fits_3, function(x) x$coefs)
pvals_3 <- sapply(fits_3, function(x) x$pvals)
resids_3 <- sapply(fits_3, function(x) x$resids)

#compare to latent factors where focal variables "unprotected"
pca2 <- prcomp(ye)
latent_vals2 <- pca2$x[,PCs_to_use]

fits_4 <- lapply(1:m, function(i){
  fit <- lm(ye[,i] ~ x[,1:k] + latent_vals2)
  return(list(
    coefs = summary(fit)$coefficients[2:(k+1),1],
    pvals = summary(fit)$coefficients[2:(k+1),4],
    resids = fit$resid
  ))
})
coefs_4 <- sapply(fits_4, function(x) x$coefs)
pvals_4 <- sapply(fits_4, function(x) x$pvals)
resids_4 <- sapply(fits_4, function(x) x$resids)

par(mfrow = c(2,4))

#pick some random outcome variable in 1:m
i = 1

#coefficients estimates are identical, perhaps as expected
plot(coefs_1[,i], coefs_2[,i], main = paste0("all variables"), 
     xlab = "coef estimates w/out latent factors", pch = 19, col = adjustcolor(1,0.2),
     ylab = "coef estimates w/ latent factors",
     xlim = range(coefs_1), ylim = range(coefs_2))
for(i in 2:m){points(coefs_1[,i], coefs_2[,i], pch = 19, col = adjustcolor(1,0.2))}
abline(0,1,col=2,lty=2,lwd=2)

#pvals for the latent factor thing are lower
plot(log10(pvals_1[,i]), log10(pvals_2[,i]), main = paste0("all variables"),
     xlab = "log10 pvals w/out latent factors", pch = 19, col = adjustcolor(1,0.2),
     ylab = "log10 pvals w/ latent factors",
     xlim = range(log10(pvals_1)), ylim =range(log10(pvals_2))); 
for(i in 2:m){points(log10(pvals_1[,i]), log10(pvals_2[,i]), pch = 19, col = adjustcolor(1,0.2))}
abline(0,1,col=2,lty=2,lwd=2)

#compare estimates to true values
plot(coefs_1[,i], b[1:k, i], main = paste0("all variables"), 
     xlab = "coef estimates from left", pch = 19, col = adjustcolor(1,0.2),
     ylab = "true values of coefficients",
     xlim = range(c(coefs_1, coefs_4)), ylim = range(b[1:k,]));
for(i in 2:m){points(coefs_1[,i], b[1:k,i], pch = 19, col = adjustcolor(1,0.2))}
abline(0,1,col=2,lty=2,lwd=2)
abline(lm(c(b[1:k,]) ~ c(coefs_1)), lwd = 2, col = 3)

#compare unprotected estimates to true values
plot(coefs_4[,i], b[1:k, i], main = paste0("all variables"), 
     xlab = "coef estimates from unprotected latent factors model", pch = 19, col = adjustcolor(1,0.2),
     ylab = "true values of coefficients",
     xlim = range(c(coefs_1, coefs_4)), ylim = range(b[1:k,])); 
for(i in 2:m){points(coefs_4[,i], b[1:k,i], pch = 19, col = adjustcolor(1,0.2))}
abline(0,1,col=2,lty=2,lwd=2)
abline(lm(c(b[1:k,]) ~ c(coefs_4)), lwd = 2, col = 3)

#poorly calibrated type-1 error 
hist(pvals_1, main = "plain \'ol observed variables")
hist(pvals_2, main = "\"latent factors\" included")
hist(pvals_3, main = "\"random noise\" included")
hist(pvals_4, main = "\"latent unprotected factors\" included")
