
#no outliers in generative model
n <- 500
nrep <- 1E4
nboots <- 1E3

pvals <- unlist(parallel::mclapply(1:nrep, function(j){
  x <- rnorm(n)
  boots <- replicate(nboots, mean(sample(x, replace = T)), T)
  mean(boots < 0)  
}, mc.cores = 12))

hist(pvals)

#outliers allowed
mixture_prop <- 0.001
sds <- c(1,50)

pvals <- unlist(parallel::mclapply(1:nrep, function(j){
  x <- rnorm(n) * sds[rbinom(n, size = 1, prob = mixture_prop) + 1]
  boots <- replicate(nboots, mean(sample(x, replace = T)), T)
  mean(boots < 0)  
}, mc.cores = 12))

hist(pvals)
