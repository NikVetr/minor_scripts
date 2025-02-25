library(corpcor)
ifelse2 <- function(bool, opt1, opt2){if(bool){return(opt1)}else{return(opt2)}}

check_psd <- function(R) all(eigen(R)$values > -1E-6)
gradleg <- function(loc, cols, labels, pos = 4, main = "", rasterize = F, border_in = NA, border_out = T, direc = "v", ...){
  
  #get metadata
  n <- length(cols)
  
  #get rect positions
  if(direc == "v"){
    locs <- data.frame(x1 = rep(loc[1], n),
                       x2 = rep(loc[2], n),
                       y1 = seq(loc[3], loc[4], length.out = n+1)[-(n+1)],
                       y2 = seq(loc[3], loc[4], length.out = n+1)[-1]
    )  
  } else {
    locs <- data.frame(x1 = seq(loc[1], loc[2], length.out = n+1)[-(length(n)+1)],
                       x2 = seq(loc[1], loc[2], length.out = n+1)[-1],
                       y1 = rep(loc[3], n),
                       y2 = rep(loc[4], n)
    )
  }
  
  #draw rects
  rect(locs$x1, locs$y1, locs$x2, locs$y2,
       col = cols, border = border_in)
  if(border_out){
    rect(loc[1], loc[3], loc[2], loc[4])
  }
  
  #draw text
  nl <- length(labels)
  if(nl != n){
    if(direc == "v"){
      locs <- data.frame(x1 = rep(loc[1], nl),
                         x2 = rep(loc[2], nl),
                         y1 = seq(loc[3], loc[4], length.out = nl+1)[-(nl+1)],
                         y2 = seq(loc[3], loc[4], length.out = nl+1)[-1]
      )  
    } else {
      locs <- data.frame(x1 = seq(loc[1], loc[2], length.out = nl+1)[-(length(nl)+1)],
                         x2 = seq(loc[1], loc[2], length.out = nl+1)[-1],
                         y1 = rep(loc[3], nl),
                         y2 = rep(loc[4], nl)
      )
    }
  }
  if(direc == "v") text(x = ifelse2(pos == 4, locs$x2, locs$x1) + 0.03, y = (locs$y1 + locs$y2) / 2, pos = pos, labels = labels)
  if(direc == "h") text(y = ifelse2(pos == 1, locs$x1, locs$x2), x = (locs$x1 + locs$x2) / 2, pos = pos, labels = labels)
  
  #draw title
  text(x = mean(loc[1:2]), y = loc[4], labels = main, pos = 3)
  
}

cmat_3 <- function(a, b, c){
  R <- diag(3)
  R[1,2] <- R[2,1] <- a
  R[1,3] <- R[3,1] <- b
  R[2,3] <- R[3,2] <- c
  R
}

try_bs <- function(a, b, c, n = 0, max_n = 20, tol = 1E-3){
  works <- check_psd(cmat_3(a, b, c))
  if(n <= max_n){
    new_n <- n + !works
    new_b <- b - 1/(2^n) * ifelse(works, 1, -1/2)
    if(abs(b-new_b) < tol){
      new_n <- max_n + 1
    }
    best_b <- try_bs(a, new_b, c, n = new_n, max_n = max_n)
  } else{
    best_b <- b - 1/(2^(n-1)) * ifelse(works, 1, -1)
  }
  return(best_b)
}
try_bs(0.5, 0.4, 0.6)
eigen(cmat_3(0.5, try_bs(0.5, 0.4, 0.6), 0.6))

try_bs(0.5, 0.4, 0.6)
0.5*0.6 - sqrt(1-0.5^2)*sqrt(1-0.6^2)

a <- 0.5
b <- 0.2
eigen(cmat_3(a, a*b + sqrt(1-a^2)*sqrt(1-b^2), b))
eigen(cmat_3(a, a*b - sqrt(1-a^2)*sqrt(1-b^2), b))

#whoops lol no numerical approx necessary!

#### analytic version ####
as <- cs <- setNames(1:10/10, 1:10/10)
bs <- -10:10/10

as <- cs <- setNames(0:200/200, 0:200/200)
da <- dc <- as.numeric(diff(as)[1])
bs <- -100:100/100

min_bs <- parallel::mclapply(as, function(ai){
  sapply(cs, function(ci){
    ai*ci - sqrt(1-ai^2)*sqrt(1-ci^2)
  })  
}, mc.cores = 12)

min_bs <- do.call(rbind, min_bs)
cols <- viridis::viridis(length(bs))
par(mar = c(5,5,2,4))
plot(as, cs, col = "white",
     xlab = latex2exp::TeX("correlation between X$_1$ & X$_2$"),
     ylab = latex2exp::TeX("correlation between X$_2$ & X$_3$"))
for(ai in seq_along(as)){
  for(ci in seq_along(cs)){
    rect(xleft = as[ai] - da/2,
         xright = as[ai] + da/2,
         ybottom = cs[ci] - dc/2,
         ytop =  cs[ci] + dc/2,
         col = cols[which.min(abs(min_bs[ai, ci] - bs))],
         border = NA)
  }
}
par(xpd = NA)
gradleg(loc = c(1.125, 1.075, 0.2, 1), labels = -5:5/5, cols = cols, main = latex2exp::TeX("r$_{1,3}$"))
text(x = 1.10, y = 1.1, labels = "min.")

#### blooming onion ####
# does the bloooming onion retrieve the same bounds?
choleskalator <- function(mat, ind){
  dim <- dim(mat)[1]
  if(ind == dim){
    return(mat)
  } else {
    ind_removed <- mgcv::choldrop(mat, ind)
    Rn_target <- c((t(mat) %*% mat[,ind])[-ind],1)  
    URn_target_sle <- backsolve(ind_removed, (Rn_target[1:(dim-1)]), transpose = T) 
    URn_target_sle <- c(URn_target_sle, sqrt(1 - crossprod(URn_target_sle))) 
    return(cbind(rbind(ind_removed, rep(0,dim-1)), URn_target_sle))
  }
}

blooming_onion_tune_chol <- function (upper_cholesky_factor, ind, varN = 0.1, betaWindow = NA, returnInOrder = T) {
  redrawBeta <- ifelse(!is.na(betaWindow), T, F)
  d <- dim(upper_cholesky_factor)[1]
  m <- d-1
  R <- matrix(0, d, d)
  chol_cor <- choleskalator(upper_cholesky_factor, ind)
  orig_order <- 1:d
  new_order <- c(orig_order[-ind], ind)
  
  R[1:m, 1:m] <- chol_cor[1:m,1:m]
  target_chol <- chol_cor[1:d,d]
  target_y <- 1 - (target_chol[d]^2)
  target_z <- target_chol[1:m] / sqrt(target_y)
  if(!is.na(betaWindow)){
    alpha <- 1 
    if(betaWindow >= 1){
      y <- rbeta(1, m/2, alpha)
    } else {
      betaBounds <- c(target_y + betaWindow, target_y - betaWindow)
      betaBounds[betaBounds > 1] <- 1
      betaBounds[betaBounds < -1] <- -1
      betaBounds_prob <- pbeta(q = betaBounds, shape1 = m/2, shape2 = alpha)
      unif_draw <- runif(1)
      beta_subprob <- betaBounds_prob[1] - betaBounds_prob[2]
      betaProb <- betaBounds_prob[2] + beta_subprob * unif_draw
      y <- qbeta(p = betaProb, shape1 = (d-1)/2, shape2 = alpha)
    }
  } else {
    y <- target_y 
  }
  z <- rnorm(m, target_z, sqrt(varN))
  z <- z/sqrt(crossprod(z)[1])
  R[1:m, m + 1] <- sqrt(y) * z
  R[m + 1, m + 1] <- sqrt(1 - y)
  
  if(redrawBeta){
    if(is.na(betaWindow)){
      
      log_prop_ratio <- 0
      
    } else {
      
      betaBounds_rev <- c(y + betaWindow, y - betaWindow)
      betaBounds_rev[betaBounds_rev > 1] <- 1
      betaBounds_rev[betaBounds_rev < -1] <- -1
      
      betaBounds_prob_rev <- pbeta(q = betaBounds_rev, shape1 = (d-1)/2, shape2 = alpha)
      beta_subprob_rev <- betaBounds_prob_rev[1] - betaBounds_prob_rev[2]
      
      log_prop_ratio <- log(abs(1 / beta_subprob_rev)) - log(abs(1 / beta_subprob)) 
      
      # does the symmetry argument for the forward and backward unit hypersphere sampling hold? not sure...
      # n_integs <- 1E2
      # range_of_integration <- seq(from = 0, to = sqrt(varN)*5, length.out = n_integs)[-1]
      # length_of_integration <- diff(range_of_integration[1:2])
      # # log_prop_ratio_normals_for <- log(sum(sapply(range_of_integration, function(multiplier) prod(dnorm(z*multiplier, target_z, sqrt(varN), log = F))))*length_of_integration)
      # ## log_prop_ratio_normals_for <- log(sum(dnorm(x = as.vector(t(t(t(range_of_integration)) %*% t(z))), mean =  rep(target_z, n_integs-1), sd =  sqrt(varN), log = F))*length_of_integration)
      # # log_prop_ratio_normals_rev <- log(sum(sapply(range_of_integration, function(multiplier) prod(dnorm(target_z*multiplier, z, sqrt(varN), log = F))))*length_of_integration)
      # ## log_prop_ratio_normals_rev <- log(sum(dnorm(x = as.vector(t(t(t(range_of_integration)) %*% t(target_z))), mean =  rep(z, n_integs-1), sd =  sqrt(varN), log = F))*length_of_integration)
      # log_prop_ratio_normals <- log_prop_ratio_normals_rev - log_prop_ratio_normals_for
      # log_prop_ratio <- log_prop_ratio + log_prop_ratio_normals
      # yes it does!! woot woot
      
    }
  } else {
    log_prop_ratio <- 0
  }
  
  if(returnInOrder){
    for(ind_rev in rep(ind, (d-ind))){
      R <- choleskalator(R, ind_rev)
    }
  }
  
  return(list(sample = R, log_prop_ratio = log_prop_ratio))
}


r <- 0.5
p <- 3
R0 <- diag(p) + r - diag(p) * r
UR0 <- chol(R0)

rlkj <- function (K, eta = 1, R0 = NULL) {
  
  if(is.null(R0)){
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
    K0 <- dim(R0)[1]
    alpha <- eta + (K - 2) / 2 - 0.5 * (K0-1)
    R <- matrix(0, K, K)
    R[1:K0, 1:K0] <- chol(R0)
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

K0 <- 5
K <- 9
R0 <- rlkj(K = K0, eta = 1)
R0 <- diag(K0) + r - diag(K0) * r
R0_mod <- R0
R0_mod[K0,K0-1] <- R0_mod[K0-1,K0] <- 0.9

R1_samp <- do.call(rbind, replicate(1E3, expr = {
  nR <-  rlkj(K = K, eta = 1, R0 = R0)
  return(nR[1:K0,K])
}, simplify = F))
apply(R1_samp, 2, mean)
hist(R1_samp[,1], breaks = 50)
hist(R1_samp[,2], breaks = 50)
hist(R1_samp[,3], breaks = 50)

#test function to modify cholesky factor
ur_digit_adjust <- function(target_r, R0 = NULL, uR0 = NULL, 
                            return_uR0_mod = F, 
                            return_new_r = T,
                            return_new_R = F,
                            return_new_uR = F){
  
  if(is.null(uR0)){
    if(is.null(R0)){
      stop("need to supply uR0 or R0")
    } else {
      uR0 <- chol(R0)
    }
  }
  uR0_mod <- uR0
  
  x <- uR0_mod[,K0-1]
  y <- uR0_mod[,K0]
  
  #constraints on bottom two digits in rightmost column or UR0
  # sum(x * y) #equals target corr
  # sum(y^2) #equals 1
  
  #need to use quadratic formula
  # for ax^2 + bx + c = 0, x = (-b {+/-} sqrt(b^2 - 4ac)) / 2a
  
  #rhs we are shooting for
  sr <- 1 - sum(uR0_mod[1:(K0-2),K0]^2)
  xyr <- target_r - sum(x[1:(K0-2)] * y[1:(K0-2)])
  
  #do algebra
  c1 <- 2 * xyr * x[K0]
  c2 <- (x[K0]^2 + x[K0-1]^2)
  c3 <- xyr^2 - sr * x[K0-1]^2
  yn <- sqrt(c3/-c2) #sign of these does not matter?
  ynm1 <- sqrt(sr - yn^2)
  uR0_mod[(K0-1):K0, K0] <- c(ynm1, yn)
  sum(uR0_mod[,K0]^2) 
  new_r <- t(uR0_mod[,K0-1]) %*% uR0_mod[,K0] #hit target!
  
  #create return object
  ifelse2 <- function(test, yes, no){if(test){return(yes)}else{return(no)}}
  out <- c(
    ifelse2(return_uR0_mod, list(mod = c(yn = yn, ynm1 = ynm1)), NULL),
    ifelse2(return_new_r, list(new_r = new_r), NULL),
    ifelse2(return_new_R, list(new_R = t(uR0_mod) %*% uR0_mod), NULL),
    ifelse2(return_new_uR, list(new_uR = uR0_mod), NULL)
  )
  
  if(length(out) == 1){
    out <- out[[1]]
  }
  
  return(out)

}

#####

r <- 0.3
p <- 5
R0 <- diag(p) + r - diag(p) * r
R0 <- rlkj(p)
r_vals <- -100:100/100
adjustments <- lapply(r_vals, function(target_r) 
  ur_digit_adjust(target_r = target_r, R0 = R0, 
                  return_new_r = T,
                  return_new_R = T,
                  return_uR0_mod = T, 
                  return_new_uR = T)
)
valid <- sapply(adjustments, function(x) !is.na(x$new_r))
adjustments <- adjustments[valid]
r_vals <- r_vals[valid]
retrieved_r_vals <- sapply(adjustments, function(x) x$new_r)
retrieved_mods <- data.frame(do.call(rbind, lapply(adjustments, function(x) x$mod)))
retrieved_R_vals <- lapply(adjustments, function(x) x$new_R)
table(sapply(retrieved_R_vals, function(x) x[K0, K0]))

#plot retrieve coef output
plot(r_vals, retrieved_r_vals, type = "l")
abline(0, 1, lwd = 10, col = adjustcolor(2, 0.3))
bounds <- round(range(retrieved_r_vals[!is.na(retrieved_r_vals)]), 4)
text(x = bounds, labels = bounds, 
     y = par("usr")[4], pos = 3, xpd = NA, col = "darkblue", font = 2)

pairs(cbind(r_vals, retrieved_mods, retrieved_r_vals))

#####

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
  } else {
    if(all(sort(indices) != c(K-1, K))){
      target_indices <- c(K-1, K)
      if(!is.null(uR)){
        uR <- choleskalator(uR, indices[2])
        uR <- choleskalator(uR, indices[1])
      } else {
        R <- reorder_R(R, indices, target_indices)
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
  x.min[K-1:0] <- c(0, sqrt_sr)
  x.max[K-1:0] <- c(sqrt_sr, 0)
  r_bounds <- c(r.min = sum(x.min * uR[,K-1]), r.max = sum(x.max * uR[,K-1]))
  return(r_bounds)
  
}

#test
p <- 5
R <- rlkj(p)
coef_bounds(R)

#lmao this only works for p=5??

#compare to numerical approach
p <- 11
R <- rlkj(p)
r_vals <- -1000:1000/1000
adjustments <- lapply(r_vals, function(target_r) 
  ur_digit_adjust(target_r = target_r, R0 = R, 
                  return_new_r = T,
                  return_new_R = T,
                  return_uR0_mod = T, 
                  return_new_uR = T)
)
valid <- sapply(adjustments, function(x) !is.na(x$new_r))
adjustments <- adjustments[valid]
r_vals <- r_vals[valid]
retrieved_r_vals <- sapply(adjustments, function(x) x$new_r)
retrieved_mods <- data.frame(do.call(rbind, lapply(adjustments, function(x) x$mod)))
retrieved_R_vals <- lapply(adjustments, function(x) x$new_R)

#plot retrieve coef output
plot(r_vals, retrieved_r_vals, type = "l")
abline(0, 1, lwd = 10, col = adjustcolor(2, 0.3))
bounds <- round(range(retrieved_r_vals[!is.na(retrieved_r_vals)]), 4)
text(x = bounds, labels = bounds, 
     y = par("usr")[4], pos = 3, xpd = NA, col = "darkblue", font = 2)
coef_bounds(R = R)
