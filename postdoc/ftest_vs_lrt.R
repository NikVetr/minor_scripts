# set.seed(1)

#functions
seq2 <- function(r, length.out, ...) seq(r[1], r[2], length.out = length.out)

#high level params
n <- 1E3 #number of individuals (eg humans) in the sample
m <- 1E4 #number of analytes / outcomes
p <- 5E0 #number of variables that systematically structure outcomes
b_sd <- 0
x_sd <- 1
R2 <- 0.1 #hmm doesn't quite hit the mark... think its an exp of var vs var of exp thing
e_var <- b_sd^2 * x_sd^2 * p * (1-R2) / R2
e_var <- 1
b <- matrix(rnorm(m*p, sd = b_sd), p, m) #true coefficients, effect of variables on each outcome
x <- matrix(rnorm(n*p), n, p) #values for each variable for each individual

#add a polynomial in there
x_foo <- x
x <- x^2
b <- b^0

#now simulate analytes
y <- x %*% b #y observed without iid noise
e <- matrix(rnorm(n*m, sd = sqrt(e_var)), n, m) #residual error for y, with ~1/3rd variance indep
ye <- y + e #the y we observe, with residual noise
x <- x_foo

fast <- T
fits <- do.call(rbind, parallel::mclapply(1:m, function(i){
  
  if(fast){
    fit <- summary(lm(ye[,i] ~ 1 + x))
    return(pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail=FALSE))
  } else {
    fit_1 <- lm(ye[,i] ~ 1 + x)
    fit_2 <- lm(ye[,i] ~ 1)
    fit_3 <- lm(ye[,i] ~ 1)
    
    c(ftest = anova(fit_1, fit_2, test = "F")$'Pr(>F)'[2],
      lrt = anova(fit_1, fit_2, test = "LRT")$'Pr(>Chi)'[2],
      wald = aod::wald.test(Sigma = vcov(fit_1), b = coef(fit_1), Terms = 1 + 1:p)$result$chi2["P"])  
  }

}, mc.cores = 12))

if(!fast){
  par(mfrow = c(3,1))
  plot(log10(fits[,1:2]))
  abline(0,1)
  
  hist(fits[,"ftest"])
  hist(fits[,"lrt"])
} else {
  par(mfrow = c(1,1))
  hist(fits)
}


#### try a new x each time ####
n <- 10000
p <- 5
nrep <- 1E3
fits <- data.frame(do.call(rbind, parallel::mclapply(1:nrep, function(i){
  
  x <- matrix(rnorm(n*p), n, p)
  y <- apply(x^2, 1, sum) + rnorm(n) * 0
  
  cov(cbind(y, x))
  
  fit <- summary(lm(y ~ 1 + x))
  return(c(ftest = pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail=FALSE),
           x = fit$coefficients[-1,4]))
  
}, mc.cores = 12)))

hist(fits$ftest.value)
alpha = 0.05
mean(fits$ftest.value < alpha)


#just a single obs
n <- 1E5
x <- rnorm(n) + 0
x <- rexp(n) * sample(c(1,-1), n, replace = T)
x <- rexp(n) * sample(c(1,-1), n, replace = T)
x <- rbinom(n, 1, 0.1)
x <- x * runif(n, 0, 9) + (1-x) * runif(n, -1, 0)
cor(x, x^2)
summary(lm(x^2 ~ x))
