n <- 10
x <- runif(n)
# x <- rep(3, n)

#pythagorean means
gmean <- function(x) exp(mean(log(x)))
hmean <- function(x) 1/(mean(1/x))
amean <- function(x) mean(x)

#other means
igmean <- function(x) log(mean(exp(x))) #-> log-sum-exp - log(n)
qmean <- function(x) sqrt(mean(x^2))

#power means
#as p -> inf, get maximum; as p -> 0, get geometric mean; as p -> 1, arithmetic; as p -> -1, harmonic
#as p -> 2, quadratic mean
pmean <- function(x, p = 2) mean(x^p)^(1/p) 

kmean <- function(x, f, i_f) i_f(mean(f(x))) #kolmogorov mean
imean <- function(x) ifelse(length(x) != 2, 
                            stop("x needs to be of length 2"), 
                            1 / exp(1) * (x[2]^x[2] / x[1]^x[1])^(1/(x[2]-x[1]))) 
bmean <- function(x) NULL #implement Bajraktarevic mean

#also frechet mean, Stolarsky mean, Chisini mean, etc.

#test them out
gmean(x)
hmean(x)
amean(x)
igmean(x)
qmean(x)
pmean(x, 1/2)

#as p -> inf, get maximum; as p -> -inf, get minimum
pmean(x, 1E2)
max(x)
pmean(x, -1E2)
min(x)

kmean(x, function(x)1/x, function(x)1/x)
kmean(x, log, exp)
kmean(x, function(x)x**2, sqrt)

#are there any means outside the range of numbers provided as input?
#do all hyperoperator means just become the min and max functions?