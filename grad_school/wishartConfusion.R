library(MCMCpack)

## I WISH TO RESOLVE WISHART CONFUSION ##

#set seed -- though you don't need to do this, since the point seems to stand across many replicates
seed(123)

#scalar degrees of freedom parameter
df <- 50 

#scale matrix
sig <- matrix(c(2,-1,-1,4),2,2)

#single random draw from Wishart with scale matrix sig and degrees of freedom df
rwish(v = df, S = sig)

#many such draws
draws <- sapply(1:2e4, function(x) rwish(v = df, S = sig))

#graphically "show" that expected matrix drawn from a Wishart so specified is the product sig*df
par(mfrow = c(2,3))
plot(density(draws[1,]), main = "Element [1,1]; Variance in Trait 1"); abline(v = sig[1,1]*df, col = 2)
plot(density(draws[2,]), main = "Element [1,2]; Covariance Between Traits 1 and 2"); abline(v = sig[1,2]*df, col = 2)
plot(density(draws[4,]), main = "Element [1,2]; Variance in Trait 2"); abline(v = sig[2,2]*df, col = 2)
plot(draws[1,], draws[2,], col=rgb(0, 0, 0, 0.15), main = "Variance in Trait 1 vs Covariance of Traits 1 and 2")
plot(draws[1,], draws[4,], col=rgb(0, 0, 0, 0.15), main = "Variance in Trait 1 vs Variance of Trait 2")
plot(draws[2,], draws[4,], col=rgb(0, 0, 0, 0.15), main = "Variance in Trait 2 vs Covariance of Traits 1 and 2")

#and yet apparently a matrix like sig*df*.9 has a larger probability density! (or less negative log-likelihood, in the case where this is a generative model)
#what gives?? 
log(dwish(sig*df, df, sig))
log(dwish(sig*df*.9, df, sig))

#also, and here's the $100 question -- sampling from the Wishart with something like rwish() returns me covariance matrices centered on df*sig
#I can trivially decompose these to either a correlation matrix (with basic linear algebra, or with a function like cov2cor), or rescale the draw
#back to the scale of the original scale matrix by dividing through by df
#as such, the Wishart seems like it would make a wonderful proposal distribution for correlation or covariance matrices
#it even comes with an intuitive tuning parameter to get that sweet -.234 acceptance probability, or whatever
#but for it to make proposals via met-hastings I need to be able to compute the proposal ratio, i.e. the ratio of the pdfs
#I can swallow that sampling covariance matrices could work perfectly well here, actually, if the above pdf confusion is resolved
#because it seems scaling by a constant shouldn't do anything to the pdf, since there'd still be a 1:1 mapping of scaled to unscaled draws
#but decomposing to a correlation matrix loses that 1:1 mapping, so a given correlation matrix can be "sampled" via the wishart in a manner that 
#comes from a draw with high probability density, or a draw with very low probability density (since the variances of the draw can be anything)
#and so knowing which density to use is tricky

riwish(v = df, S = sig)




## note, also this seems consistent with what the internet says, e.g. https://www.statlect.com/probability-distributions/wishart-distribution#hid3

# no need to look below this, just sketched it out for my own satisfaction
#construct a covariance matric -- in this case a correlation matrix -- with dimensionality D
D <- 5; 
L_Sig <- matrix(0,D,D);
L_Sig[1,1] <- 1;
for (i in 2:D) {
  bound <- 1;
  for (j in 1:(i-1)) {
    L_Sig[i,j] <- runif(1, -sqrt(bound), sqrt(bound));
    bound <- bound - L_Sig[i,j]^2;
  }
  L_Sig[i,i] <- sqrt(bound);
}
Sig <- L_Sig %*% t(L_Sig);


Sig
n <- 10000
matrices <- array(dim = c(D,D,n))
for(i in 1:n){
  a <- rmvnorm(1, sigma = Sig) 
  matrices[,,i] <- t(a) %*% a
}

#the mean matrix is very close to sig

meanMatrix <- function(matrices){
  rows <- dim(matrices)[1]
  cols <- dim(matrices)[2]
  samples <-  dim(matrices)[3]
  meanMatrix <- matrix(0, nrow = rows, ncol = cols)
  for ( i in 1:rows ){
    for( j in i:cols ){
      meanMatrix[i,j] <- mean(matrices[i,j,])
    }
  }
  meanMatrix <- meanMatrix + t(meanMatrix) - diag(diag(meanMatrix))
  return(meanMatrix)
}

Sig
meanMatrix(matrices)
