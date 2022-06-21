#check understanding
n <- 1E5
nexp <- 1E2
hist(rgamma(n = n, shape = nexp, rate = 1))
hist(replicate(n, sum(rexp(nexp, 1))), add = T, col = adjustcolor(2, 0.2))

#read in data
d <- readLines("~/Desktop/SA-2021-predictions.txt", warn = F)
d <- do.call(rbind, strsplit(d, ": "))
probs <- as.integer(gsub("%", "", d[,2]) ) / 100
happened <- c(8, 14, 17, 18, 21, 27:30, 38, 40, 45, 46, 48:50, 54, 60, 61, 65:67, 69, 75, 79, 80, 83, 86, 92, 94, 95, 100, 103:105, 108)

#evaluate log loss
sum_log_loss <- sum(log(probs[happened])) + sum(log(1-probs[-happened]))
pgamma(-sum_log_loss, shape = length(probs), rate = 1)
pgamma(-length(probs)*log(0.5), shape = length(probs), rate = 1)
  

mode_conc <- function(w, k){
  if(k > 2){
    plot(0:1000/1000, dbeta(0:1000/1000, w*(k-2)+1, (1-w)*(k-2) + 1), type = "l")  
  } else {
    plot(0:1000/1000, dbeta(0:1000/1000, w*k, (1-w)*k), type = "l")
  }
  abline(v = w)
}


mode_conc(0.1, 1.9999)

#get some 3x3 correlation matrix
corx <- 0.7 + diag(3) * (1-0.7)

#define rotation matrices in 3D
rotx <- function(t) matrix(rbind(c(1,0,0), c(0, cos(t), -sin(t)), c(0,sin(t),cos(t))), 3, 3)
roty <- function(t) matrix(rbind(c(cos(t),0,sin(t)), c(0, 1, 0), c(cos(t),-sin(t),0)), 3, 3)
rotz <- function(t) matrix(rbind(c(cos(t),-sin(t),0), c(sin(t), cos(t), 0), c(0,0,1)), 3, 3)

#get final rotation matrix
rad <- runif(3, 0, 2*pi)
rotmat <- rotx(rad[1]) %*% rotz(rad[2]) %*% rotz(rad[3])

#get eigendecomp of correlation matrix
eigen <- eigen(corx)
v <- eigen$vectors
l <- eigen$values

#get orig correlation matrix back 
v %*% diag(l) %*% t(v)

#fiddle with the eigenvalues?
newl <- sort(runif(3, 0, 1), T)
newl <- newl / sum(newl) * 3
v %*% diag(newl) %*% t(v) #nope :[

#fiddle with the eigenvectors?
newv <- rotmat %*% v

#still orthonormal?
newv %*% t(newv)
apply(newv^2, 2, sum)

#does it work? 
newv %*% diag(l) %*% t(newv) #nah

#gamma distributions
mu <- 0.11
s2 <- 1E3
gamma_params <- function(mu, s2) c(mu^2/s2, mu / s2)
ab <- gamma_params(mu,s2)
mean(rgamma(1E4, ab[1], ab[2]))
var(rgamma(1E4, ab[1], ab[2]))
