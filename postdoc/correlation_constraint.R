library(corpcor)
ifelse2 <- function(bool, opt1, opt2){if(bool){return(opt1)}else{return(opt2)}}
source("~/scripts/minor_scripts/postdoc/3-way_correlation_constraint.R")

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
  
  K0 <- dim(R0)[1]
  
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
    print(paste0("using indices: (", indices[1], ", ", indices[2], ")"))
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

#test
p <- 5
R <- rlkj(p)
coef_bounds(R)

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

#looks to work?
range(retrieved_r_vals)
coef_bounds(R = R)

#explicit test
nR <- R
nR[p,p-1] <- nR[p-1,p] <- coef_bounds(R = R)[2] + 1E-9
eigen(nR)$value


nR[p,p-1] <- nR[p-1,p] <- coef_bounds(R = R)[1] - 1E-9
eigen(nR)$value

#wait this does not work... or rather it works for the upper and not the lower bound
#probably to do with the sqrt, +/- thing
#can maybe assess by multiplying uR by -1?


#### conditional correlations ####

#now let's figure out the conditional distribution
frac_end <- function(p){
  z <- rnorm(p)
  z <- z/sqrt(crossprod(z)[1])
  z[p-1]^2 / sum(z[(p-1):p]^2)
}

#ratio is F distributed
hist(replicate(1E4, frac_end(5)))

#### peeling the onion ####

#initial params
p <- 11
R <- rlkj(p)
R0 <- R[1:(p-1), 1:(p-1)]
eta <- 1

#forward direction sampling
p0 <- p-1
alpha <- eta + (p - 2) / 2 - 0.5 * (p0-1)
uR <- matrix(0, p, p)
uR[1:p0, 1:p0] <- chol(R0)
p2a <- p - p0
if (p2a > 0) {
  for (m in p0:(p0 + p2a - 1)) {
    alpha <- alpha - 0.5
    y <- rbeta(1, m/2, alpha)
    z <- rnorm(m, 0, 1)
    z <- z/sqrt(crossprod(z)[1]) #unit hypersphere rescaling
    uR[1:m, m + 1] <- sqrt(y) * z
    uR[m + 1, m + 1] <- sqrt(1 - y)
  }
}
crossprod(uR)

#unit hypersphere test
sample_hypersphere <- function(m){
  z <- rnorm(m, 0, 1)
  z/sqrt(crossprod(z)[1])
}
p <- 10
n <- 1E5
breaks = -100:100/100
hist(replicate(n, sample_hypersphere(p)[p]), freq = F, ylim = c(0,2), breaks = breaks)
hist(sqrt(rbeta(n, shape1 = 1/2, shape2 = (p-1)/2)) * (rbinom(n, 1, 0.5)*2-1), 
     add = T, freq = F, col = adjustcolor(2,0.5), breaks = breaks)
curve(2 * abs(x) * dbeta(x^2, shape1 = 1/2, shape2 = (p-1)/2) / 2, 
      from = -1, to = 1, n = 512, add = TRUE, col = 3, lwd = 3)

#go backwards from uR
uR <- chol(R)
x <- uR[,p]
pcol <- uR[,p]
m <- p - 1
pcol[1:(m-1)] #these are fixed by values at other coefficients
sr <- 1 - sum(pcol[1:(m-1)]^2) #final two coefficients have sum of squares remainder sr

#try getting y first? and can then check if z_m has appropriate distribution
nsamp <- 1E4
alpha <- eta + (p - 2) / 2 - 0.5 * m
rbeta(1, m/2, alpha) #is the untruncated distribution for y
#but if uR_{p,p} is sqrt(1 - y)
#then sqrt(1 - y)^2 cannot be more than sr
#in other words sqrt(1 - y)^2 < sr
#so 1 - y < sr and 1 - sr < y
#or rather y > 1 - sr
y_lb <- 1-sr
y_q_lb <- qbeta(p = y_lb, shape1 = m / 2, shape2 = alpha)
y_runif <- runif(nsamp)
y_q <- y_q_lb + (1-y_q_lb) * y_runif
y_samp <- pbeta(q = y_q, shape1 = m / 2, shape2 = alpha)

#wait, x_p is not likely to be all of sr is it? the z_m will not have 
#same expectation of z_{1,2,...,m-1}?
#although maybe that is just the constraint imposed by R / uR?

x_p <- sqrt(1 - y_samp)
x_m <- sqrt(sr - x_p^2)

# x.new[p] <- x_p
# x.new[m] <- x_m
# 
# # r.old <- sum(x * uR[,m]) #same as R[m,p]
# r.new <- sum(x.new * uR[,m])
r.new <- sum(x[1:(m-1)] * uR[1:(m-1),m]) + x_m * uR[m,m]

hist(r.new)
abline(col = 2, lwd = 3, v = R[m,p])

range(r.new)
coef_bounds(R)

#let's compare to a numerical approx?
# nsamp_rlkj <- 1E5
# R0 <- R[1:m, 1:m]
# uR0 <- chol(R0)
# last_corrs <- do.call(cbind, lapply(1:nsamp_rlkj, function(i) 
#   rlkj(K = p, eta = eta, uR0 = uR0)[1:p,p]
# ))
# 
# distances_for_other_corrs <- apply(last_corrs[1:(m-1),] - R[1:(m-1),p], 2, function(x) sum(x^2))
# nearby_corrs <- which(distances_for_other_corrs < quantile(distances_for_other_corrs, 0.0001))
# r.new.numeric <- last_corrs[m,nearby_corrs]
# hist(r.new, col = 2, breaks =  -20:20/20, freq = F)
# hist(r.new.numeric, breaks = -20:20/20, add = T, col = adjustcolor(1,0.5), freq = F)
# 
# plot(last_corrs[1:(m-1),nearby_corrs[1]],
# R[1:(m-1),p])

#ok, I think this is not working because the conditional of y and z is jointly distributed
#aka Pr(y,z|x_{1:m-1}) has some complicated interdependence
#let's try sampling with Stan instead

#### check unit hypersphere rotation ####


check_rot_hyper <- function(k){
  x <- rnorm(k)
  
  #sample directly from unit hypersphere
  xk <- x / sqrt(crossprod(x)[1])
  
  #retrieve infrasphere sample
  xkm1 <- x[-k] / sqrt(crossprod(x[-k])[1])
  
  #resample kth element via rotation
  theta <- runif(1, min = 0, max = pi/2)
  y <- c(sin(theta) * xkm1, cos(theta))
  
  #via normal rescaling
  z <- c(xkm1 * sqrt(rchisq(1, k-1)), rnorm(1))
  z <- z / sqrt(crossprod(z)[1])
  
  #via marginal distribution
  w_beta <- rbeta(1, 1/2, (k-1)/2)
  w <- c(xkm1 * (1-w_beta), sqrt(w_beta) * sample(c(-1,1), 1))
  
  return(c(from_rot = y[k], from_chi = z[k], from_beta = w[k], direct = xk[k], implied_theta = acos(xk[k])))
}

k <- 10
last_elements <- data.frame(do.call(rbind, replicate(1E4, check_rot_hyper(k), simplify = F)))
breaks = -20:20/20
hist(last_elements$direct, breaks = breaks)
hist(last_elements$from_rot, breaks = breaks, col = adjustcolor(2, 0.5), add = T)
hist(last_elements$from_chi, breaks = breaks, col = adjustcolor(3, 0.5), add = T)
hist(last_elements$from_beta, breaks = breaks, col = adjustcolor(4, 0.5), add = T)
hist(last_elements$implied_theta) #range is always 0, pi
mean(last_elements$implied_theta) #mean is always pi/2
sd(last_elements$implied_theta) #sd varies with k
#when k = 50, sd = 1/7
#when k = 5, sd = 1/2
#when k = 100, sd = 1/sqrt(99)?
#when k = 100, sd = 1/sqrt(99)?
#when k = 4, sd = ...almost 1/sqrt(3)
var(last_elements$implied_theta) 

#is this a stretched beta?
var(last_elements$implied_theta / pi) * 2
var(rbeta(1E4, k/2, k/2))

#hmm maybe not

#maybe we can just sample correlation in Stan directly... and hardcode the constraint / bounds?

#### Stan conditional sampling ####
library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

#sample full dimension correlation matrix
p <- 10
R <- rlkj(p)
uR <- chol(R)

#build Stan model
stan_loc <- "~/scripts/minor_scripts/postdoc/conditional_unit-hypersphere-x-beta.stan"
stan_loc_direct <- "~/scripts/minor_scripts/postdoc/conditional_corrmat_direct.stan"
stan_program <- paste0(readLines(stan_loc), collapse = "\n")
stan_program_direct <- paste0(readLines(stan_loc_direct), collapse = "\n")

do_Stan_sampling <- F
if(do_Stan_sampling){
  
mod <- cmdstan_model(stan_loc)
mod_direct <- cmdstan_model(stan_loc_direct)

#build data object
dat <- list(m = p - 1,
            eta = 1,
            x_obs = uR[1:(p-2),p],
            rpm1 = uR[1:(p-1),p-1])

dat_direct <- list(m = p,
            eta = 1,
            R_obs = R,
            r_bounds = unlist(coef_bounds(R)))


fit <- mod$sample(chains = 4, iter_sampling = 5E3, iter_warmup = 5E3, data = dat, 
                  parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, 
                  refresh = 100, init = 0.1)
fit_direct <- mod_direct$sample(chains = 4, iter_sampling = 5E3, iter_warmup = 5E3, data = dat_direct, 
                  parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, 
                  refresh = 100, init = 0.1)
# summ <- fit$summary()
samps <- as.data.frame(as_draws_df(fit$draws()))
samps_direct <- as.data.frame(as_draws_df(fit_direct$draws()))
hist(samps$r, breaks = 100)
hist(samps_direct$r, breaks = 100)
}

#### range frac comparison ####
smoooth_EVs <- function(R, n){
  evs <- eigen(R)
  p <- dim(R)[1]
  logEVs <- log(evs$values[1:(n-1)])
  spline_fit <- smooth.spline(1:(n-1), logEVs)
  logEVs <- c(logEVs, predict(spline_fit, n:p)$y)
  newEVs <- exp(logEVs)
  newEVs <- newEVs / sum(newEVs) * p
  newR <- (evs$vectors) %*% diag(newEVs) %*% t(evs$vectors)
  return(cov2cor(newR))
}


estimate_conditional <- F
p <- 100
R <- rlkj(p)
r <- 0.8
R <- diag(p) + r - diag(p) * r
uR <- chol(R)
n <- 20
x <- matrix(rnorm(n*p), ncol = p, nrow = n) %*% uR
sample.R.true <- cor(x)
sample.R <- as.matrix(Matrix::nearPD(cor(x), corr = T)$mat)
I_weight <- 0.0
sample.R <- sample.R * (1-I_weight) + diag(p) * I_weight
# sample.R <- smoooth_EVs(cor(x), n = n)
sample.uR <- chol(sample.R)
all.indices <- t(combn(1:p, 2))
r.bounds <- data.frame(do.call(rbind, lapply(1:choose(p,2), function(i){
  
  #try approximating the conditional distribution (if p < n, or else constraint is exact)
  indices <- all.indices[i,]
  # print(paste0(i, ": ", paste0(indices, collapse = ", ")))
  target_indices <- c(p-1, p)
  uR_shuff <- choleskalator(choleskalator(sample.uR, indices[2]), indices[1])
  r.obs <- sample.R[indices[1], indices[2]]
  
  #build data object
  conditional_quantile <- NA
  if(estimate_conditional){
    dat <- list(m = p - 1,
                eta = 1,
                x_obs = uR_shuff[1:(p-2),p],
                rpm1 = uR_shuff[1:(p-1),p-1])
    
    #fit model
    sink(tempfile())
    fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat, 
                      parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, 
                      refresh = 100, init = 0.1)
    sink()
    samps <- as.data.frame(as_draws_df(fit$draws()))
    
    #assess sample value in this distribution
    conditional_quantile <- mean(r.obs > samps$r)  
  }
  
  
  #now put it all together
  out <- c(row = indices[1],
    col = indices[2],
    r.true = R[indices[1], indices[2]], 
    r.obs.est = r.obs,
    r.obs.samp = sample.R.true[indices[1], indices[2]],
    coef_bounds(uR = sample.uR, indices = indices),
    cond_q = conditional_quantile
  )
  return(out)
})))

#get a few extra useful variables
r.bounds$r.obs <- r.bounds$r.obs.samp
r.bounds$r.prop <- (r.bounds$r.obs - r.bounds$r.min) / (r.bounds$r.max - r.bounds$r.min)
r.bounds$sample.error <- r.bounds$r.obs - r.bounds$r.true
r.bounds$sample.error.z <- atanh(r.bounds$r.obs) - atanh(r.bounds$r.true)
r.bounds$deviation.from.mid <- r.bounds$r.obs - (r.bounds$r.max + r.bounds$r.min) / 2
r.bounds$deviation.from.mid.z <- atanh(r.bounds$r.obs) - atanh((r.bounds$r.max + r.bounds$r.min) / 2)

#neat, these are the same
plot(atanh(r.bounds$r.obs), atanh((r.bounds$r.max + r.bounds$r.min) / 2))

#do some plotting
logit <- function(p) log(p/(1-p))
# plot(abs(r.bounds$r.prop - 0.5) * 2, abs(r.bounds$sample.error))
# cor(abs(r.bounds$r.prop - 0.5), abs(r.bounds$sample.error))
plot(r.bounds$r.prop, r.bounds$sample.error,
     main = paste0("pearson's r = ", 
                   round(cor(r.bounds$r.prop, r.bounds$sample.error), 2)),
     pch = 19, col = adjustcolor(1,0.2),
     xlab = "relative location in range of PSD constrained interval",
     ylab = "error of sample correlation from true correlation"
)
abline(h=0,col=2,lty=2,lwd=2, xpd=F)
abline(v=0.5,col=2,lty=2,lwd=2, xpd=F)

if(estimate_conditional){
  plot(r.bounds$cond_q, r.bounds$sample.error,
       main = paste0("pearson's r = ",
                     round(cor(r.bounds$r.prop, r.bounds$sample.error), 2)),
       pch = 19, col = adjustcolor(1,0.2),
       xlab = "relative location in range of PSD constrained interval",
       ylab = "error of sample correlation from true correlation"
  )
  abline(h=0,col=2,lty=2,lwd=2, xpd=T)
  abline(v=0.5,col=2,lty=2,lwd=2, xpd=T)
}


#### iteratively updating observed matrix? ####

#generate data
p <- 8
R <- rlkj(p)
# r <- 0.5
# R <- diag(p) + r - diag(p) * r
uR <- chol(R)
n <- 5
x <- matrix(rnorm(n*p), ncol = p, nrow = n) %*% uR
sample.R <- as.matrix(Matrix::nearPD(cor(x), corr = T)$mat)
# sample.R <- smoooth_EVs(cor(x), n = n)

#estimate coef bounds
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


bounded_sample.R <- all_coef_bounds(sample.R)
bounded_sample.R$wb <- bounded_sample.R$ub - bounded_sample.R$lb
bounded_sample.R$pb <- (sample.R - bounded_sample.R$lb) / bounded_sample.R$wb
bounded_sample.R$dir <- 0.5 - bounded_sample.R$pb 
cor(R[upper.tri(R)], (sample.R + bounded_sample.R$dir / 10)[upper.tri(R)])
cor(R[upper.tri(R)], sample.R[upper.tri(R)])
