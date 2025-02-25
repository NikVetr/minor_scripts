library(corpcor)

#### functions ####

ifelse2 <- function(bool, opt1, opt2){if(bool){return(opt1)}else{return(opt2)}}
source("~/scripts/minor_scripts/postdoc/3-way_correlation_constraint.R")
source("~/repos/polylines/R/functions.R")
rlkj <- function (K, eta = 1, R0 = NULL, uR0 = NULL) {
  
  if(is.null(R0) & is.null(uR0)){
    alpha <- eta + (K - 2) / 2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K)
    R[1, 1] <- 1
    R[1, 2] <- r12
    R[2, 2] <- sqrt(1 - r12^2)
    if (K > 2) {
      for (m in 2:(K - 1)) {
        alpha <- alpha - 0.5
        y <- rbeta(1, m/2, alpha)
        z <- rnorm(m, 0, 1)
        z <- z/sqrt(crossprod(z)[1]) #unit hypersphere rescaling
        R[1:m, m + 1] <- sqrt(y) * z
        R[m + 1, m + 1] <- sqrt(1 - y)
      }
    }
    
  } else {
    if(is.null(uR0)){
      uR0 <- chol(R0)
    }
    K0 <- dim(uR0)[1]
    alpha <- eta + (K - 2) / 2 - 0.5 * (K0-1)
    R <- matrix(0, K, K)
    R[1:K0, 1:K0] <- uR0
    k2a <- K - K0
    if (k2a > 0) {
      for (m in K0:(K0 + k2a - 1)) {
        alpha <- alpha - 0.5
        y <- rbeta(1, m/2, alpha)
        z <- rnorm(m, 0, 1)
        z <- z/sqrt(crossprod(z)[1]) #unit hypersphere rescaling
        R[1:m, m + 1] <- sqrt(y) * z
        R[m + 1, m + 1] <- sqrt(1 - y)
      }
    }
  }
  
  return(crossprod(R))
  
}

choleskalator <- function(mat, ind){
  
  #given corr matrix 'R' and its UR cholesky factor 'mat'
  #cheaply rearranges 'mat' to obtain 'out', such that Q = t(out) %*% out 
  #and Q corresponds to R with the variable in location 'ind'
  #now in the last row and column
  dim <- dim(mat)[1]
  if(ind == dim){
    return(mat)
  } else {
    ind_removed <- mgcv::choldrop(mat, ind)
    Rn_target <- c((t(mat) %*% mat[,ind])[-ind],1)  
    URn_target_sle <- backsolve(ind_removed, (Rn_target[1:(dim-1)]), transpose = T) 
    URn_target_sle <- c(URn_target_sle, sqrt(1 - crossprod(URn_target_sle))) 
    out <- cbind(rbind(ind_removed, rep(0,dim-1)), 
                 URn_target_sle)
    colnames(out) <- rownames(out) <- NULL
    return(out)
  }
}

reorder_R <- function(R, indices, target_indices){
  K <- dim(R)[1]
  x_inds <- 1:K
  y_inds <- 1:K
  x_inds[c(indices[1], target_indices[1])] <- 
    x_inds[c(target_indices[1], indices[1])]
  y_inds[c(indices[2], target_indices[2])] <- 
    y_inds[c(target_indices[2], indices[2])]
  R <- R[y_inds, y_inds]
  R <- R[x_inds, x_inds]
  return(R)
}


#now let's write a function to retrieve the bounds of the coefficient
#it looks like y_{n-1} determines this -- 
#when it is 0 and y_{n} fills the remainder of the simplex, we get the lower bound
#when it fills the remainder and y_{n} is 0, r is at its maximum value

coef_bounds <- function(R = NULL, uR = NULL, indices = NULL){
  
  #retrieve dimension of matrix
  if(is.null(uR)){
    if(is.null(R)){
      stop("need to supply uR or R")
    } else {
      K <- dim(R)[1]
    }
  } else {
    K <- dim(uR)[1]
  }
  
  #rearrange matrix or ur-cholesky factor (via the choleskalator) to simplify compute
  if(is.null(indices)){
    indices <- c(K-1, K)
    # print(paste0("using indices: (", indices[1], ", ", indices[2], ")"))
  } else {
    if(!all(sort(indices) == c(K-1, K))){
      target_indices <- c(K-1, K)
      if(!is.null(uR)){
        uR <- choleskalator(uR, indices[2])
        uR <- choleskalator(uR, indices[1])
      } else {
        R <- reorder_R(R, indices = indices, target_indices = target_indices)
      }
    }
    
  }
  
  #retrieve upper right cholesky factor if we do not yet have it
  if(is.null(uR)){
    uR <- chol(R)
  }
  
  #extract last column to evaluate bounds
  x.min <- x.max <- uR[,K]
  sr <- 1-sum(uR[1:(K-2), K]^2) #simplex remainder
  sqrt_sr <- sqrt(sr)
  partial_sum <- sum(uR[1:(K-2),K-1] * uR[1:(K-2),K])
  Km1_coef <- uR[K-1,K-1]
  
  r_bounds <- c(r.min = partial_sum - abs(Km1_coef) * sqrt_sr,
                r.max = partial_sum + abs(Km1_coef) * sqrt_sr)
  
  # x.min[K-1:0] <- c(0, sqrt_sr)
  # x.max[K-1:0] <- c(sqrt_sr, 0)
  # r_bounds <- c(r.min = sum(x.min * uR[,K-1]), 
  #               r.max = sum(x.max * uR[,K-1]))
  return(r_bounds)
  
}
all_coef_bounds <- function(R){
  p <- dim(R)[1]
  uR <- chol(R)
  all.indices <- t(combn(1:p, 2))
  all.bounds <- data.frame(do.call(rbind, lapply(1:choose(p,2), function(i){
    
    indices <- all.indices[i,]
    target_indices <- c(p-1, p)
    uR_shuff <- choleskalator(choleskalator(uR, indices[2]), indices[1])
    
    out <- c(row = indices[1],
             col = indices[2],
             r.obs = R[indices[1], indices[2]], 
             coef_bounds(uR = uR, indices = indices)
    )
    return(out)
  })))
  
  lb  <- ub <- diag(p)
  lb[as.matrix(all.bounds[,c("row", "col")])] <- all.bounds$r.min
  ub[as.matrix(all.bounds[,c("row", "col")])] <- all.bounds$r.max
  lb <- lb + t(lb) - diag(p)
  ub <- ub + t(ub) - diag(p)
  return(list(lb = lb, ub = ub))
}


#test
p <- 500
rij <- 0.3
R <- diag(p) + rij - diag(p) * rij
coef_bounds(R)

#### get bound widths for variety of p and rij ####

#compound symmetry case
nps <- 50
max_p <- 200
ps <- unique(floor(exp(seq(from = 1, to = log(max_p), length.out = nps))))
ps <- ps[ps>2]
nps <- length(ps)
nrij <- 50
eps <- 1E-5
mcprint <- function(...){system(sprintf('printf "%s"', paste0(..., collapse="")))}
constraint_widths_cs <- parallel::mclapply(ps, function(p){
  mcprint(paste0(p, " "))
  rijs <- seq(from = -1/(p-1) + eps, to = 1-eps, length.out = nrij)
  bound_widths <- unlist(lapply(rijs, function(rij){
    R <- diag(p) + rij - diag(p) * rij
    bounds <- coef_bounds(R)
    return(diff(bounds))
  }))
  return(data.frame(p = p, rij = rijs, width = bound_widths))
}, mc.cores = 8)

#sample variation case
nps <- 50
max_p <- 200
ps <- unique(floor(exp(seq(from = 1, to = log(max_p), length.out = nps))))
ps <- ps[ps>2]
nps <- length(ps)
nrij <- 50
eps <- 1E-5
nreps <- 1
nx <- max_p + 10
mcprint <- function(...){system(sprintf('printf "%s"', paste0(..., collapse="")))}
constraint_widths_samp <- lapply(ps, function(p){
  mcprint(paste0(p, " "))
  rijs <- seq(from = -1/(p-1) + eps, to = 1-eps, length.out = nrij)
  bound_widths <- unlist(parallel::mclapply(rijs, function(rij){
    R <- diag(p) + rij - diag(p) * rij
    L <- t(chol(R))
    dbounds <- sapply(1:nreps, function(ri) {
      x <- matrix(rnorm(p*nx), nx, p)
      Rs <- cor(t(L %*% t(x)))
      bounds <- all_coef_bounds(Rs)
      return(mean((bounds$ub - bounds$lb)[upper.tri(Rs)]))  
    })
    return(mean(dbounds))
  }, mc.cores = 14))
  return(data.frame(p = p, rij = rijs, width = bound_widths))
})

#### plotting ####
cols <- viridisLite::viridis(nps)
plot(constraint_widths_cs[[1]]$rij, constraint_widths_cs[[1]]$width, 
     type = "l", xlim = c(-0.5,1), ylim = c(0,2), col = cols[1],
     xlab = "off-diagonal correlation (compound symmetry)",
     ylab = "width of marginal PSD-preserving interval", frame.plot = F)
for(i in 2:nps){
  lines(constraint_widths_cs[[i]]$rij, constraint_widths_cs[[i]]$width, 
        type = "l", xlim = c(-1,1), ylim = c(0,2), col = cols[i])
}
gradleg(loc = c(1,1.1,1,2), cols = cols, 
        labels = ps[ceiling(seq(1, nps, length.out = 10))], 
        rasterize = T, main = "Dim.", border_out = T)

plot(constraint_widths_samp[[1]]$rij, constraint_widths_samp[[1]]$width, 
     type = "l", xlim = c(-0.5,1), ylim = c(0,2), col = cols[1],
     xlab = "average off-diagonal correlation (100 replicates)",
     ylab = "width of marginal PSD-preserving interval", frame.plot = F)
for(i in 2:nps){
  lines(constraint_widths_samp[[i]]$rij, constraint_widths_samp[[i]]$width, 
        type = "l", xlim = c(-1,1), ylim = c(0,2), col = cols[i])
}

gradleg(loc = c(1,1.1,1,2), cols = cols, 
        labels = ps[ceiling(seq(1, nps, length.out = 10))], 
        rasterize = T, main = "Dim.", border_out = T)


#### test identity sample constraint ####
p <- 200
nx <- 200
rij <- 1E-5
R <- diag(p) + rij - diag(p) * rij
L <- t(chol(R))
x <- matrix(rnorm(p*nx), nx, p)
Rs <- cor(t(L %*% t(x)))
bounds <- coef_bounds(Rs)

