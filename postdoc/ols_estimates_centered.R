pval <- function(fit){
  summary(fit)$coefficients[,4]
}

foo <- function(n, b2_fac = 1){
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  e1 <- rnorm(n)
  e2 <- rnorm(n)
  
  b <- runif(1)
  b <- rt(1, df = 2)
  b <- rnorm(1)
  b <- rexp(1)
  b <- abs(rnorm(1))
  b <- qnorm(runif(1, 0.4, 0.6)) #truncated normal
  b <- rnorm(1, sd = 0.1) #small normal
  
  a <- rnorm(1)
  
  fits <- list(
    lm((a + b * x1 + e1)~x1), 
    lm((a + b * x1 + e2)~x1),
    lm((a + b * x2 + e1)~x2),
    lm((a + b2_fac * b * x2 + e2)~x2)
    )
  
  out <- c(b = b, 
           est. = do.call(rbind, lapply(fits, coef))[,2],
           pval. = do.call(rbind, lapply(fits, pval))[,2],
           pval.b = 0,
           est.b = b)
  return(out)
}

dev.off()

nrep <- 2E4
b2_fac <- 1.5
out <- do.call(rbind, parallel::mclapply(1:nrep, function(i) foo(100, b2_fac = b2_fac), mc.cores = 16))
pairs(out[1:min(1E3, nrep),1:5], pch = 19, col = adjustcolor(1, 0.5))

#ols regression coefficients
cout <- cov(out[,1:5])
cout %*% diag(1/diag(cout))

#hmm, unbiased estimates do *NOT* imply they come off 1-to-1 line, cos if you have 2 normal dists next to each other the line going through them is ~0
#depends on SE of estimate

tls_estimate <- function(x, y){
  v <- prcomp(cbind(x,y))$rotation
  fit <- v[2,1]/v[1,1]
  fit <- c(mean(y) - fit * mean(x), fit)
  fit
}

#how does tls perform?
tls_slopes <- sapply(setNames(1:5, colnames(out)[1:5]), function(yi){
  sapply(setNames(1:5, colnames(out)[1:5]), function(xi){
    fit <- tls_estimate(out[,xi],out[,yi])
    fit[2] #2nd element is the slope
  })
})

tls_slopes
(tls_slopes + t(tls_slopes)) / 2
max(abs(t(tls_slopes) - 1 / tls_slopes)) #also they are each other's inverse

tls_slopes[c("est.1","est.4"),c("est.1","est.4")]

#tls slope guaranteed when b is normal (exponential? student's t?) and both sides equally powered, 
#OR when b is anything (finite variance RV?) and both sides are infinitely powered?
#ah no, it is just at the limit of infinite power? or just under equal power
#expected value of tls slope under equal power = 1, pretty sure

#now find index 2 & 5's sampling dist (independent x and e)
nbootstrap <- 5E3
# bstrap_slopes <- t(replicate(n = nbootstrap, {
#   inds <- sample(1:nrow(out), replace = T)
#   tls_estimate(out[inds,2],out[inds,5])
#   }))[,2]
bstrap_slopes_eqpow <- do.call(rbind, parallel::mclapply(1:nbootstrap, function(i){
  inds <- sample(1:nrow(out), replace = T)
  tls_estimate(out[inds,"est.1"],out[inds,"est.4"])
}, mc.cores = 16))[,2]
quantile(bstrap_slopes_eqpow, probs = c(0.025, 0.975))

bstrap_slopes_uneqpow <- do.call(rbind, parallel::mclapply(1:nbootstrap, function(i){
  inds <- sample(1:nrow(out), replace = T)
  tls_estimate(out[inds,"b"],out[inds,"est.1"])
}, mc.cores = 16))[,2]
quantile(bstrap_slopes_uneqpow, probs = c(0.025, 0.975))

#now see what filtering by pvalues does to these estimates
alpha = 0.05
comp.inds <- c("b",4)
passing <- which(p.adjust(out[,paste0("pval.", comp.inds[1])], method = "BH") < alpha &
  p.adjust(out[,paste0("pval.", comp.inds[2])], method = "BH") < alpha)
sigout <- out[passing, c(paste0("est.", comp.inds[1]), paste0("est.", comp.inds[2]))]

cov(sigout) %*% diag(1/diag(cov(sigout)))
plot(sigout)
tls_estimate(sigout[,paste0("est.", comp.inds[1])], sigout[,paste0("est.", comp.inds[2])])
bstrap_slopes_sigfilt <- do.call(rbind, parallel::mclapply(1:nbootstrap, function(i){
  inds <- sample(1:nrow(sigout), replace = T)
  tls_estimate(sigout[inds,1],sigout[inds,2])
}, mc.cores = 16))[,2]
quantile(bstrap_slopes_sigfilt, probs = c(0.025, 0.975))

#is tls estimator very stable? nope@
# nx <- 1E2
# x <- cbind(rnorm(nx), rnorm(nx))
# hist(do.call(rbind, parallel::mclapply(1:nbootstrap, function(i){
#   inds <- sample(1:nrow(x), replace = T)
#   tls_estimate(x[inds,1], x[inds,2])
# }, mc.cores = 16))[,2], breaks = 10000, xlim = c(-10,10))
# 
# hist(do.call(rbind, parallel::mclapply(1:nbootstrap, function(i){
#   inds <- sample(1:nrow(x), replace = T)
#   coef(lm(x[inds,1] ~ x[inds,2]))
# }, mc.cores = 16))[,2], breaks = 100, xlim = c(-0.1,0.1))

#when b2_fac =/= 0, very weakly powered to detect deviation in double significant hits
#but decently powered to detect variation overall
#but in both cases, very strongly powered to detect variation in power between groups!