e <- 2.71828182845
iterMCMC <- 1e6
sampInt <- 200
up <- 0
down <- 0
same <- 0
curr <- 35
samples <- vector(length = iterMCMC/sampInt)
testLLs <- vector(length = iterMCMC)

postMean <- 0
postSD <- 100

currLL <- log(dnorm(x = curr, mean = postMean, sd = postSD))

for (i in 1:iterMCMC){
  if(i%%10000==0){print(c(i, currLL, curr))}
  # test <- curr + rnorm(n = 1, mean = 0, sd = 10)
  test <- curr * rnorm(n = 1, mean = 1, sd = .01) #BIASED TO SHRINK
  testLL <- log(dnorm(x = test, mean = postMean, sd = postSD))
  testLLs[i] <- testLL
  if(testLL > currLL){
    currLL <- testLL
    curr <- test
    up <- up + 1
  } else if (runif(n = 1, min = 0, max = 1) < e^(testLL - currLL)){
    currLL <- testLL
    curr <- test
    down <- down + 1
    # print("lower")    }
  } else {same <- same + 1}
  if(i %% sampInt == 0){
    samples[[i/sampInt]] <- curr
  }
}
samples <- samples[-(1:((iterMCMC/sampInt)/4))]
dens(samples)
mean(samples); sd(samples)


## BM TR?

CovMatrix <- function(num){ #generate a random covariance matrix with size of num x num
  sigmas <-  (runif(num, min = 3, max = 5.5))
  Rho <- matrix(runif(num^2, min=-.2, max=.4), nrow=num)
  diag(Rho) <- rep(1,num)
  Rho[lower.tri(Rho)] = t(Rho)[lower.tri(Rho)]
  return((diag(sigmas) %*% Rho %*% diag(sigmas)))
}
Cov <- CovMatrix(10)
state1 <- rnorm(10, 5, 3)
state2 <- rnorm(10, 5, 3)
dmvnorm(x = state1, mean = state2, sigma = Cov)
dmvnorm(x = state2, mean = state1, sigma = Cov)
