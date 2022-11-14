#### simple model ####

#simulation specifications
n_fit <- 2E3
n_test <- 1E3

#top level parameters
p <- 5
mean_pvarx <- 0.1
conc_pvarx <- 2

#correlation between residuals of x
rij_resid <- 0.1
R <- diag(p) * (1-rij_resid) + rij_resid

#proportion variance of X mappable to Y
random_betas_variances <- F
if(random_betas_variances){
  pvarx <- rbeta(p, mean_pvarx * conc_pvarx, (1-mean_pvarx)*conc_pvarx)
  B <- rnorm(p)
  B <- B / mean(abs(B))
} else {
  B <- rep(1, p)
  pvarx <- rep(mean_pvarx, p)
}

#simulate training set
y_fit <- rnorm(n_fit)
x_fit_resid <-  matrix(rnorm(n_fit * p), n_fit, p) %*% chol(R)
x_fit <- sapply(1:p, function(i){B[i] * y_fit * sqrt(pvarx[i] * B[i]^2) + x_fit_resid[,i] * sqrt((1-pvarx[i]) * B[i]^2)})

#simulate test set
y_test <- rnorm(n_test)
x_test_resid <- matrix(rnorm(n_test * p), n_test, p) %*% chol(R)
x_test <- sapply(1:p, function(i){B[i] * y_test * sqrt(pvarx[i] * B[i]^2) + x_test_resid[,i] * sqrt((1-pvarx[i]) * B[i]^2)})

#fit model and compute predictions on test set
fit <- lm(y_fit ~ x_fit)
preds <- data.frame(y_test = y_test, y_pred = cbind(1, x_test) %*% t(t(coefficients(fit))))
(var(preds$y_test) - var(preds$y_test - preds$y_pred)) / var(preds$y_test)
cor(preds)^2
plot(preds)

#### quick test -- are resids covariances described by Schur complement? yep ####
n <- 1E3
p <- 10
rij <- 0.1
ind0 <- 3
inds <- ind0:p
R <- diag(p) * (1-rij) + rij

x <- matrix(rnorm(n * p), n, p)
# x <- matrix(runif(n * p), n, p) %*% chol(R)
# x <- matrix(rbinom(n * p, 1, 0.2), n, p) %*% chol(R)
# 
# myCop <- copula::mvdc(copula=copula::normalCopula(param=R[upper.tri(R)], dim = p, dispstr = "un"), 
#               margins=rep("binom", p),
#               paramMargins=lapply(1:p, function(x) list(size=2, prob=0.25)))
# x <- copula::rMvdc(n, myCop)

x <- x %*% diag(1 / sqrt(apply(x, 2, var))) %*% solve(chol(cor(x))) %*% chol(R)
cov(lm(x[,1] ~ x[,inds])$resid, x[,-c(1,inds)])

covx <- cov(x) #this makes everything within sample, not out of sample -- use true covar to get out of sample
A <- covx[-inds, -inds]
C <- covx[inds, inds]
B <- covx[-inds, inds]
cond_covx <- (A - B %*% solve(C) %*% t(B))
cond_covx[1,-c(1,inds)]

cond_covx[1,1]
var(lm(x[,1] ~ x[,inds])$resid)

summary(lm(x[,1] ~ x[,inds]))$r.squared
1 - cond_covx[1,1] / covx[1,1]

#the covariances of the Y|X are given by the schur complement of 'Y' block (A) in the overall covariance matrix
#which also describes the covariances of the residuals of a linear model with unobserved covariates that covary with predictors 
#would be wishart distributed if cov is a sample covariance

#### alternate model ####

#simulation specifications
n_fit <- 2E3
n_test <- 1E3

#top level parameters
p <- 100
mean_pvarx <- 0.1
conc_pvarx <- 2
U_pvarxresid <- 0.0

#proportion variance of X mappable to Y
random_betas_variances <- F
if(random_betas_variances){
  pvarx <- rbeta(p, mean_pvarx * conc_pvarx, (1-mean_pvarx)*conc_pvarx)
  B <- rnorm(p)
  B <- B / mean(abs(B))
} else {
  B <- rep(1, p)
  pvarx <- rep(mean_pvarx, p)
}

#simulate training set
y_fit <- rnorm(n_fit)
U_fit <- rnorm(n_fit)
x_fit_resid <-  matrix(rnorm(n_fit * p), n_fit, p) * sqrt(1-U_pvarxresid) + U_fit * sqrt(U_pvarxresid)
x_fit <- sapply(1:p, function(i){B[i] * y_fit * sqrt(pvarx[i] * B[i]^2) + x_fit_resid[,i] * sqrt((1-pvarx[i]) * B[i]^2)})

#simulate test set
y_test <- rnorm(n_test)
U_test <- rnorm(n_test)
x_test_resid <- matrix(rnorm(n_test * p), n_test, p) * sqrt(1-U_pvarxresid) + U_test * sqrt(U_pvarxresid)
x_test <- sapply(1:p, function(i){B[i] * y_test * sqrt(pvarx[i] * B[i]^2) + x_test_resid[,i] * sqrt((1-pvarx[i]) * B[i]^2)})

#fit model and compute predictions on test set
fit <- lm(y_fit ~ x_fit)
y_pred <- cbind(1, x_test) %*% t(t(coefficients(fit)))
preds <- data.frame(y_test = y_test, y_pred = y_pred, y_resid = y_test - y_pred)
cor(preds)^2
plot(preds$y_pred, preds$y_test, pch=19,col=adjustcolor(1, 0.5)); abline(0,1,lty=2,lwd=2,col=2)

(var(preds$y_test) - var(preds$y_resid)) / var(preds$y_test)
1-exp(diff(log(diag(cov(preds[,c("y_test", "y_resid")])))))

cov(U_test, preds$y_resid)
cov(U_test, preds$y_test)

#### compute conditional covariance matrix explicitly (assuming random_betas_variances = F) ####
# n <- 1E4; a <- rnorm(n); b <- a + rnorm(n); summary(lm(a ~ b))
p <- 100
mean_pvarx <- 0.1
conc_pvarx <- 2
U_pvarxresid <- 0.0

total_cov <- diag(p+2)
colnames(total_cov) <- rownames(total_cov) <- c("y", "U", paste0("x_", 1:p))
covx <- diag(p) * (1- (mean_pvarx + (1-mean_pvarx) * U_pvarxresid)) + (mean_pvarx + (1-mean_pvarx) * U_pvarxresid)
covyx <- rep(mean_pvarx, p)
covUx <- rep((1-mean_pvarx) * U_pvarxresid, p)
total_cov[-c(1:2), -c(1:2)] <- covx
total_cov[1, -c(1:2)] <- total_cov[-c(1:2), 1] <- covyx
total_cov[2, -c(1:2)] <- total_cov[-c(1:2), 2] <- covUx

inds <- 3:(p+2)
A <- total_cov[-inds, -inds]
C <- total_cov[inds, inds]
B <- total_cov[-inds, inds]
cond_cov <- (A - B %*% solve(C) %*% t(B))
1 - cond_cov[1,1] / total_cov[1,1]



#### another quick test -- confounding model ####
n <- 1E3
npred <- 10
p <- npred + 1
rij_x1 <- 0.2
var_prop <- 0 #a number between 0 and 1 representing the distance between max and min variance in rij_x1
rij_xn1e <- 0.0
rij_xn1e <- min(rij_xn1e, 1 - rij_x1^2)
ind0 <- 2
inds <- ind0:p
R <- diag(p)
R[1,inds] <- R[inds,1] <- rij_x1 #corr of outcome with predictors
R[inds,inds] <- diag(p-ind0+1)*(1-rij_x1^2-rij_xn1e) + rij_x1^2 + rij_xn1e #corr of predictors due to above relation, plus extra

x <- matrix(rnorm(n * p), n, p)
x <- x %*% diag(1 / sqrt(apply(x, 2, var))) %*% solve(chol(cor(x))) %*% chol(R)
cov(lm(x[,1] ~ x[,inds])$resid, x[,-c(1,inds)])

covx <- cov(x) #this makes everything within sample, not out of sample -- use true covar to get out of sample
A <- covx[-inds, -inds]
C <- covx[inds, inds]
B <- covx[-inds, inds]
if(is.null(dim(A))){A <- t(A)}
if(is.null(dim(B))){
  if(length(inds) == 1){
    B <- t(t(B))
  } else {
    B <- t(B)
  }
  
}
if(is.null(dim(C))){C <- t(C)}

cond_covx <- (A - B %*% solve(C) %*% t(B))
cond_covx[1,-c(inds)]

cond_covx[1,1]
var(lm(x[,1] ~ x[,inds])$resid)

summary(lm(x[,1] ~ x[,inds]))$r.squared
1 - cond_covx[1,1] / covx[1,1]


#### put it in a function and plot ####
diag2 <- function(x){
  if(length(x) > 1){
    return(diag(x))
  } else {
    if(abs(x%%1)>1E-6){
      return(as.matrix(x))
    } else{
      return(diag(x))
    }
  }
}

R2_from_Cov <- function(p, y_rij, u_rij_prop, long_way = F){
  
  #adjust supplied u_rij if impossible
  if(p == 0){return(0)}
  u_rij <- u_rij_prop * (1 - y_rij^2)
  
  if(long_way){
    #generate covariance matrix
    k <- p + 2
    inds <- 3:k
    R <- diag(k)
    R[1,inds] <- R[inds,1] <- y_rij
    R[2,inds] <- R[inds,2] <- u_rij
    R[inds,inds] <- diag(p)*(1-y_rij^2-u_rij) + y_rij^2 + u_rij  
    
    A <- R[-inds, -inds]
    B <- R[-inds, inds]
    C <- R[inds, inds]
    Ci <- solve(C) #slow method for general inverses, vs rank 1 update
  } else {
    A <- diag(2)
    B <- rbind(rep(y_rij, p), rep(u_rij, p))
    Ci <- fastmatrix::sherman.morrison(a = diag2(rep(1 / (1 - y_rij^2 - u_rij), p)),
                                       b = rep(y_rij^2 + u_rij, p),
                                       d = rep(1, p), inverted = T)
  }
  
  #solve schur complement
  cond_R <- A - B %*% Ci %*% t(B)
  
  #return R2
  1 - cond_R[1,1] / R[1,1]
}

np <- sort(unique(c(0, round(sqrt(2)^(0:18)))))
y_rijs <- 1:9/10
u_rijs <- 0:99/100

R2s <- lapply(u_rijs, function(ur){
  sapply(y_rijs, function(yr){
    sapply(np, function(pi){
      R2_from_Cov(pi, yr, ur)
      })
    })
})
R2s <- abind::abind(R2s, along = 3, new.names = list(p = np, y_r = y_rijs,  u_r = u_rijs))

y_rij <- y_rijs[9]
plot(u_rijs, R2s["0",as.character(y_rij),], type = "l", ylim = c(0,1))
for(pi in np){
  lines(u_rijs, R2s[as.character(pi),as.character(y_rij),])
}


#### test relationships in 3-variable system ####
n <- 2E6
p <- 2
y <- rnorm(n)
u <- rnorm(n)
r_yx <- 0.85
r_ux_prop_remain <- 0.1
r_xx <- r_yx^2 + r_ux_prop_remain * (1 - r_yx^2)
b1 <- 1
b2 <- 0
b3 <- sqrt((1 / r_yx)^2 - 1)
#three way system
b1 <- 1
b2 <- sqrt(r_xx / r_yx^2 - 1)
b3 <- sqrt(1 / r_yx^2 - 1 - b2^2)
r_ux <- b2 / sqrt(b1^2 + b2^2 + b3^2)
r_ux <- sqrt(r_xx - r_yx^2) #in terms of just the r_ijs
x <- replicate(p, y * b1 + u * b2 + rnorm(n) * b3)
x <- x %*% diag(1 / apply(x, 2, sd))
cor(cbind(y, u, x))
c(r_yx, r_ux, r_xx)
# cov(cbind(y, u, x))
#corrs among the ’x’s are sqrt(r_y,x_i)
#u’s variance is constrained to everything residual of y’s

# R2_from_Cov(p, r_yx, r_ux_prop_remain)
cor(lm(y ~ x)$resid, lm(u ~ x)$resid) #partial correlation
cor(lm(y ~ x)$resid, u) #semi-partial correlation

#get cond cor
k <- p + 2
inds <- 3:k
R <- diag(k)
R[1,inds] <- R[inds,1] <- r_yx
R[2,inds] <- R[inds,2] <- r_ux 
R[inds,inds] <- diag(p) * (1 - r_xx) + r_xx

A <- R[-inds, -inds]
B <- R[-inds, inds]
C <- R[inds, inds]
Ci <- solve(C)
cond_R <- A - B %*% Ci %*% t(B)

cond_R[1,2] / sqrt(cond_R[1,1]) / sqrt(cond_R[2,2])
cor(lm(y ~ x)$resid, lm(u ~ x)$resid) #partial correlation

cond_R[1,2] / sqrt(cond_R[1,1])
cor(lm(y ~ x)$resid, u) #semi-partial correlation
