library(matrixcalc) # for matrixcalc::is.positive.semi.definite

simulateLKJMatrix = function(eta, dim) {

  # initialize the matrices
  P = matrix(0, dim, dim)
  S = diag(dim)

  # initialize beta
  beta = eta + (dim - 1) / 2

  for(k in 1:(dim - 1)) {

    # decrement beta
    beta = beta - 1/2

    for(i in (k+1):dim) {
      P[k,i] = rbeta(1, beta, beta) * 2 - 1
      p = P[k,i]
      if(k > 1) {
        for(l in (k-1):1) {
          p = p * sqrt((1-P[l,i]^2)*(1-P[l,k]^2)) + P[l,i] * P[l,k]
        }
      }
      S[i,k] = p
      S[k,i] = p
    }

  }

  return(S)

}

eta = 0.01
k   = 10
R   = lapply(1:10000, function(x) simulateLKJMatrix(eta, k))

ipsd = sapply(R, is.positive.semi.definite)

mean(ipsd)

