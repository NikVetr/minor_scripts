p <- 200
r <- 0.1
R <- diag(p) + r - diag(p) * r
L <- t(chol(R))
nrep <- 100

check_eval <- function(p, n, L = NULL){
  if(is.null(L)){
    r <- 0.5
    R <- diag(p) + r - diag(p) * r
    L <- t(chol(R))
  } 
  x <- t(L %*% matrix(rnorm(p*n), p, n))
  return(c(svar = var(x[,1]),
           seigen1 = eigen(cov(x))$values[1]))
}

var_vs_eigen <- t(replicate(100, check_eval(p = p, n = 100, L = L)))
var_vs_eigen <- apply(var_vs_eigen, 2, function(x) x / mean(x))
apply(var_vs_eigen, 2, var)

check_evec <- function(p, n, L = NULL){
  if(is.null(L)){
    r <- 0.5
    R <- diag(p) + r - diag(p) * r
    L <- t(chol(R))
  } 
  x <- t(L %*% matrix(rnorm(p*n), p, n))
  return(eigen(cor(x))$vector[,1])
}

r <- 0.1
R <- diag(p) + r - diag(p) * r
L <- t(chol(R))
nrep <- 100
evecs <- replicate(nrep, check_evec(p, 100, L = L))

find_angle <- function(vec1, vec2){
  rads <- acos(sum(vec1 * vec2))
  rads <- min(abs(rads), abs(pi - rads))
  return(rads / 2 / pi * 360)
}

find_rpair_angle <- function(evecs){
  inds <- sample(1:ncol(evecs), 2)
  return(find_angle(evecs[,inds[1]], evecs[,inds[2]]))
}
npairs <- 1E3
hist(replicate(npairs, find_rpair_angle(evecs)), breaks = 100)

generate_hypersphere_points <- function(dim, n_points, L = NULL) {
  points <- matrix(rnorm(n_points * dim), ncol = dim)
  if(!is.null(L)){
    points <- t(L %*% t(points))
  }
  points <- points / sqrt(rowSums(points^2))
  return(points)
}

randvecs <- t(generate_hypersphere_points(p, nrep, L = L))
hist(replicate(npairs, find_rpair_angle(randvecs)), breaks = 100)
