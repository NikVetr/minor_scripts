library(mvtnorm) 
library(pbivnorm)

rlkj <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}

condMVN <- function(means, cov, obs, inds){
  A <- as.matrix(cov[-inds, -inds])
  C <- as.matrix(cov[inds, inds])
  B <- as.matrix(cov[-inds, inds])
  if(dim(B)[1] != length(means)){B <- t(B)}
  condCov <- A - B%*%solve(C)%*%t(B)
  condMean <- as.vector(means[-inds] + B%*%solve(C)%*%(obs - means[inds]))
  return(list(means = condMean, cov = condCov))
}

testfun <- function(){

dimr <- 10
# lower = runif(dimr, -1, 0) 
lower = rep(-Inf, dimr)
upper = runif(dimr, 0, 1)
# upper = rep(Inf, dimr)
disp <- rnorm(dimr)
lower <- lower + disp
upper <- upper + disp
R <- rlkj(dimr)

true_integ <- pmvnorm(mean = rep(0,dimr), lower = lower, upper = upper, corr = R, algorithm = GenzBretz())[1]
true_integ

n_shuffles <- 10
approx <- rep(NA, n_shuffles)
for(iter in 1:n_shuffles){
  integs <- rep(NA, dimr)
  shuffled_inds <- sample(1:dimr, size = dimr, replace = F)
  R_sub <- R[shuffled_inds,shuffled_inds]
  ind <- 1
  means <- rep(0,dimr)
  while(dim(R_sub)[1] != 1){
    nsq <- 50
    qs <- seq(1E-6, pnorm(q = upper[ind], mean = means[1], sd = sqrt(R_sub[1,1])), length.out = nsq)
    locs <- qnorm(qs, mean = means[1], sd = sqrt(R_sub[1,1]))
    ps <- dnorm(x = locs, mean = means[1], sd = sqrt(R_sub[1,1]))
    integs[ind] <- pnorm(q = upper[ind], mean = means[1], sd = sqrt(R_sub[1,1]), log.p = T)
    conditionals <- lapply(locs, function(loc) condMVN(means = means, cov = R_sub, inds = 1, obs = loc))
    means <- apply(do.call(rbind, lapply(1:nsq, function(i) conditionals[[i]]$means * ps[i])), 2, sum) / sum(ps)
    R_sub <- conditionals[[1]]$cov * ps[1]
    for(i in 2:nsq){
      R_sub <- R_sub + conditionals[[i]]$cov * ps[i]
    }
    R_sub <- R_sub / sum(ps)
    ind <- ind + 1
  }
  integs[dimr] <- pnorm(q = upper[ind], mean = means, sd = sqrt(R_sub[1,1]), log.p = T)
  approx[iter] <- sum(integs)
}

return(c(approx = mean(approx, na.rm = T), true = log(true_integ)))

}

test <- replicate(5, testfun())
plot(t(test))
mod <- lm(test[2,] ~ test[1,])
abline(mod)
mod

#####


testfunc <- function(){
  
  dimr <- 10
  # lower = runif(dimr, -1, 0) 
  lower = rep(-Inf, dimr)
  upper = runif(dimr, 0, 1)
  # upper = rep(Inf, dimr)
  disp <- rnorm(dimr)
  lower <- lower + disp
  upper <- upper + disp
  R <- rlkj(dimr)
  
  true_integ <- pmvnorm(mean = rep(0,dimr), lower = lower, upper = upper, corr = R, algorithm = GenzBretz())[1]
  true_integ
  
  integs <- rep(NA, dimr)
  means <- rep(0,dimr)
  for(ind in 1:dimr){
    nsq <- 100
    runifs <- replicate(nsq, runif(n = dimr-1, 0, 1))
    qs <- diag(pnorm(q = upper[-ind], mean = means[-ind], sd = sqrt(diag(R)[-ind]))) %*% runifs 
    locs <- qnorm(qs) #assumes mean 0, sd 1
    ps <- sapply(1:nsq, function(x) dmvnorm(locs[,x]))
    draw_inds <- sample(1:nsq, replace = T, size = nsq, prob = ps)
    draw_inds <- table(draw_inds)
    conditionals <- lapply(as.integer(names(draw_inds)), function(draw_ind) 
      condMVN(means = means, cov = R, inds = (1:dimr)[-ind], obs = locs[,draw_ind]))
    mean_for_ind <- sum(sapply(1:length(draw_inds), function(i) conditionals[[i]]$means * draw_inds[i])) / nsq
    var_for_ind <- sum(sapply(1:length(draw_inds), function(i) conditionals[[i]]$cov * draw_inds[i])) / nsq
    
    integs[ind] <- pnorm(q = upper[ind], mean = var_for_ind, sd = sqrt(var_for_ind), log.p = T)
  }
  

  return(c(approx = sum(integs), true = log(true_integ)))
  
}

testfunc()

test <- replicate(50, testfunc())
plot(t(test))
mod <- lm(test[2,] ~ test[1,])
abline(mod)
mod
