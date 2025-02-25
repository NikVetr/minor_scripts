library(cmdstanr)

#### test intersect function ####

# 1) Compile the Stan model for the intersection function
mod_file <- "~/scripts/minor_scripts/postdoc/eigencmat-intersectHyperellipsoid.stan"
mod <- cmdstan_model(mod_file)

# 2) Prepare some test data
p_test <- 50
L <- sort(unitvec(p_test)^2, T)^2
Lambda_test <-   L <- L / sum(L) * p
y_test <- 0.8        # already a "Beta^(1/k)" value in [0,1]
sign_flip_test <- 1  # or -1

# Put them in a list for Stan
dat_list <- list(
  p = p_test,
  Lambda = Lambda_test,
  y = y_test,
  sign_flip = sign_flip_test
)

# 3) Run the model in fixed_param mode
fit <- mod$sample(
  data = dat_list,
  seed = 123,
  chains = 1,
  iter_warmup = 0,
  iter_sampling = 1,    # let's get 5 draws
  thin = 1,
  refresh = 0,          # suppress console output
  fixed_param = TRUE    # crucial for a model with no posterior
)

# 4) Extract the draws for x
draws_x <- fit$draws("x", format = "matrix")
print(draws_x)
sum(draws_x[1,]^2)
sum(draws_x[1,]^2 * Lambda_test)

#### test conditional intersection function ####

# 1) Create the .stan file
mod_file <- "~/scripts/minor_scripts/postdoc/eigencmat-sampleIntersection.stan"
# writeLines(...  # or just save manually

# 2) Compile the Stan model
mod <- cmdstan_model(mod_file)

# 3) Provide data
V_curr_test <- V[,1:(i-1), drop = F] #from cmat-eigenconstraints.R
B_curr_test <- B_curr
p_test <- nrow(V_curr_test)
pcurr_test <- ncol(V_curr_test)
L_test <- L
y_test <- 0.7      # a "Beta^(1/(prem-2))" in [0,1], e.g. 0.7
sign_flip_test <- 1  # or -1

dat_list <- list(
  p       = p_test,
  pcurr   = pcurr_test,
  L       = L_test,
  V_curr  = V_curr_test,
  B_curr  = B_curr_test,
  y       = y_test,
  sign_flip = sign_flip_test
)

# 4) Fit the model in fixed_param mode because there's no real posterior
fit <- mod$sample(
  data = dat_list,
  seed = 123,
  chains = 1,
  iter_warmup = 0,
  iter_sampling = 1,
  fixed_param = TRUE, 
  refresh = 0
)

# 5) Extract results
draws_out <- fit$draws("out", format="matrix")

# Each row of draws_out is one draw from the 'generated quantities' matrix out.
# The columns represent out[,1], out[,2], ..., out[,prem], flattened.
# So if prem = 4, out has 4 columns. 'out' is p x prem = 6 x 4 = 24 total elements, 
# flattened in column-major order => out[1,1], out[2,1],..., out[p,1], out[1,2],..., out[p,2], etc.

# For example, let's recover the last draw as a matrix
last_draw <- draws_out[nrow(draws_out), ]  # a 24-length numeric if prem=4
# We'll reshape into p x prem
p <- p_test
prem <- p_test - pcurr_test
out_matrix <- matrix(last_draw, nrow=p, ncol=prem)
new_vec <- out_matrix[, 1]
B_new   <- out_matrix[, 2:prem]
cat("Newly sampled vector:\n"); print(new_vec)
sum(new_vec^2)
sum(new_vec^2*L)
sapply(1:pcurr_test, function(pi) sum(new_vec * V_curr_test[,pi]))
cat("Updated basis:\n"); print(B_new)
sapply(1:ncol(B_new), function(pi) sum(new_vec * B_new[,pi]))


#### test construction function ####

#compile the model
mod_file <- "~/scripts/minor_scripts/postdoc/eigencmat-constructCorrmat.stan"
mod <- cmdstan_model(mod_file)

#simulate data
p_test <- 20
L <- sort(unitvec(p_test)^2, T)
L_test <-   L <- L / sum(L) * p_test

# We need p-1 = 4 draws for y in [0,1]
y_array_test <- runif(p_test-1, 0, 1)

# We need p-1 = 4 sign flips
sign_flips_test <- sample(c(-1,1), p_test-1, T)    # each in {-1, +1}

# better yet store as a vector
theta_values <- numeric(0)   # will append to this
for (i in 1:(p_test - 1)) {
  subsize <- (p_test - i - 1)  # = p - i - 1
  if (subsize > 0) {
    # sample a random direction in R^subsize, normalized
    v <- rnorm(subsize)
    v <- v / sqrt(sum(v^2))
    theta_values <- c(theta_values, v)
  }
}

# Put them in a list
dat_list <- list(
  p               = p_test,
  L               = L_test,
  y_array         = y_array_test,
  sign_flip_array = sign_flips_test,
  theta_values    = theta_values
)

# fit the Stan model
fit <- mod$sample(
  data = dat_list,
  seed = 123,
  chains = 1,
  iter_warmup = 0,
  iter_sampling = 1,
  fixed_param = TRUE,
  refresh = 0
)

mat_draws <- fit$draws("C_samp", format="matrix")
C_samp_last <- matrix(mat_draws[nrow(mat_draws), ], p_test, p_test)
diag(C_samp_last)
range(abs(eigen(C_samp_last)$values - L))

#### sample correlation matrices ####

mod_file <- "~/scripts/minor_scripts/postdoc/eigencmat-sampleCorrmat.stan"
mod <- cmdstan_model(mod_file)

# Supply only p
dat_list <- list(p=10)

fit <- mod$sample(
  data = dat_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000
)

# Now you can extract C_samp or any other parameter
summ <- fit$summary()
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])
samps <- fit$draws("C_samp")
samps <- fit$draws("C_samp")
# Extract dimensions
iter_sampling <- dim(samps)[1]  # Number of posterior samples per chain
chains <- dim(samps)[2]         # Number of MCMC chains
p <- sqrt(dim(samps)[3])        # Infer p (since C_samp is stored as p x p flattened)
C_array <- array(samps, dim = c(iter_sampling * chains, p, p))
C_array <- aperm(C_array, c(2, 3, 1))  # Rearrange to (p, p, total_samples)

#### fit correlation matrices ####
mod_file <- "~/scripts/minor_scripts/postdoc/eigencmat-fitCorrmat-debug.stan"
# mod_file <- "~/scripts/minor_scripts/postdoc/eigencmat-fitCorrmat-debug-fixtheta.stan"
mod <- cmdstan_model(mod_file)

# simulate data
n <- 100
p <- 5
r_ij <- 0.3
Sigma <- matrix(r_ij, p, p)
diag(Sigma) <- 1
x_data <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
Sigma_samp <- cor(x_data)
dat_list <- list(
  n = n,
  p = p,
  x = x_data,
  incl_ll = 0 #if 1 then add to target, if 0 do not (sample from prior)
)
if(grepl("fixtheta", mod_file)){
  n_angles <- ((p - 3) * (p - 2)) / 2
  theta_angles <- runif(n_angles, 0, pi)
  dat_list[length(dat_list) + 1] <- list(theta_angles)
  names(dat_list)[length(dat_list)] <- "theta_angles"
}

# fit
fit <- mod$sample(
  data = dat_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,  
  iter_warmup = 500,
  iter_sampling = 1000, 
  refresh = 10, 
  adapt_delta = 0.9, 
  max_treedepth = 12
)

#inspect
summ <- fit$summary("C_samp")
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])
samps <- fit$draws("C_samp")
# Extract dimensions
iter_sampling <- dim(samps)[1]  # Number of posterior samples per chain
chains <- dim(samps)[2]         # Number of MCMC chains
p <- sqrt(dim(samps)[3])        # Infer p (since C_samp is stored as p x p flattened)
C_array <- array(samps, dim = c(iter_sampling * chains, p, p))
cmats <- aperm(C_array, c(2, 3, 1))  # Rearrange to
par(mfrow = c(p,p), mar = c(1,1,1,1) * 1.1)
for(i in 1:p){
  for(j in 1:p){
    if(i <= j){
      plot.new()
    } else {
      hist(cmats[i,j,], breaks = -100:100/100, 
           main = paste0("(", i, ", ", j, ")"),
           xlab = "", ylab = "", yaxt = "n", xaxt = "n")
      segments(-1,0,1,0, xpd = NA)
      segments(-1,0,-1,par("usr")[4], xpd = NA)
      segments(par("usr")[2],0,par("usr")[2],par("usr")[4], 
               xpd = NA)
      abline(v=Sigma[i,j], col = 2, lwd = 2)
      abline(v=Sigma_samp[i,j], col = 3, lwd = 2)
      abline(v=0,lty=2)
    }
  }
}


#### fit fixed param model ####
mod_file <- "~/scripts/minor_scripts/postdoc/eigencmat-fitCorrmat-debug-fixparams.stan"
mod <- cmdstan_model(mod_file)

# simulate data
n <- 10
p <- 8
r_ij <- 0.0
Sigma <- matrix(r_ij, p, p)
diag(Sigma) <- 1
x_data <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
Sigma_samp <- cor(x_data)
#   simplex[p] raw_dirichlet;
#   array[p - 1] real<lower=-1, upper=1> signed_y_array;
#   vector<lower=0, upper=pi()>[((p - 3) * (p - 2)) %/% 2] theta_angles;
raw_dirichlet <- c(gtools::rdirichlet(1, rep(1, p)))
signed_y_array <- runif(p-1, -1, 1)
n_angles <- ((p - 3) * (p - 2)) / 2
theta_angles <- runif(n_angles, 0, pi)
theta_angles <- theta_angles * runif(n_angles, 0.99, 1.01)

#construct data input object
dat_list <- list(
  n = n,
  p = p,
  x = x_data,
  incl_ll = 1, #if 1 then add to target, if 0 do not (sample from prior)
  raw_dirichlet = raw_dirichlet,
  signed_y_array = signed_y_array,
  theta_angles = theta_angles
)

# fit
fit <- mod$sample(
  data = dat_list,
  seed = 123,
  chains = 1,
  parallel_chains = 1,  
  iter_warmup = 0,
  iter_sampling = 1, 
  refresh = 10, 
  fixed_param = T
)

#inspect
samps_theta <- data.frame(posterior::as_draws_df(fit$draws("theta_values")))
samps_theta <- samps_theta[,-c(ncol(samps_theta)-0:2)]
samps_C <- fit$draws("C_samp")
# Extract dimensions
iter_sampling <- dim(samps_C)[1]  # Number of posterior samples per chain
chains <- dim(samps_C)[2]         # Number of MCMC chains
p <- sqrt(dim(samps_C)[3])        # Infer p (since C_samp is stored as p x p flattened)
C_array <- array(samps_C, dim = c(iter_sampling * chains, p, p))
cmats <- aperm(C_array, c(2, 3, 1))  # Rearrange to

# foo <- cmats[,,1]
# bar <- samps_theta
cmat_diffs <- (cmats[,,1] - foo)[upper.tri(foo)]
theta_diffs <- unlist(samps_theta - bar)
par(mfrow = c(3,1))
hist(cmat_diffs, breaks = -10:10/10)
hist(theta_diffs)
abs_cmat_diffs <- (abs(cmats[,,1]) - abs(foo))[upper.tri(foo)]
plot(abs(cmat_diffs), abs(abs_cmat_diffs))

#PROBLEM -- small differences in theta can cause the sign of
#specific correlation elements to flip!
#TODO check if this is also the case in the R code