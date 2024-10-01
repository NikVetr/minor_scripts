library(cmdstanr)
library(posterior)
library(lavaan)

#### sim data ####

#specify high level parameters
standardize_vars <- F
p = 50
n = 200
prop_yvar = 0.5 #proportion non-aleatoric variance (ie, not in error term)
prop_zvar = 0.2
sd_b_xz_scale <- 0.2
sd_b_xz = sd_b_xz_scale / sqrt(p)
b_yz = 4
sd_u_scale <- 0 #relative confounding effect of unobserved_variable
sd_u <- sd_u_scale * sqrt(p)

#set stage to simulate correlated counts in {0,1,2}
rop <- 0.8
rs <- c(1, rop^(1:p / 3))
R <- outer(1:p, 1:p, FUN = function(i, j, rs) rs[abs(i - j) + 1], rs = rs)
cholR <- chol(R)
p_freqs <- rbeta(p, 2, 2)
p_liab <- qnorm(p_freqs)

#actually simulate these
x1 <- Reduce("+", lapply(1:2, function(i){
  liabs <- matrix(rnorm(p * n), nrow = n) %*% cholR
  liabs
}))

x2 <- Reduce("+", lapply(1:2, function(i){
  liabs <- matrix(rnorm(p * n), nrow = n) %*% cholR
  liabs
}))

x3 <- Reduce("+", lapply(1:2, function(i){
  liabs <- matrix(rnorm(p * n), nrow = n) %*% cholR
  liabs
}))

#further randomize effects of x?
# x1 <- x1 %*% matrix(rnorm(p^2), p, p)
# x2 <- x2 %*% matrix(rnorm(p^2), p, p)

#specify mediator, outcome, and exposure relationships
b_xy <- rnorm(p)
b_xz <- rnorm(p) * sd_b_xz #when sd_b_xz = 0, no horizontal pleiotropy 

### simulate data for ###
### z <- x -> y -> z  ###

#simulate y for pop 1
mu_y1 <- c(x1 %*% t(t(b_xy)))
u1 <- rnorm(n) * sd_u
e_y1 <- rnorm(n) * c(sqrt(var(mu_y1 + u1) / prop_yvar * (1-prop_yvar)))
y1 <- mu_y1 + u1 + e_y1

#simulate y and z for pop 2
mu_y2 <- c(x2 %*% t(t(b_xy)))
u2 <- rnorm(n) * sd_u
e_y2 <- rnorm(n) * c(sqrt(var(mu_y2 + u2) / prop_yvar * (1-prop_yvar)))
y2 <- mu_y2 + u2 + e_y2
mu_z2 <- y2 * b_yz + c(x2 %*% t(t(b_xz)))
e_z2 <- rnorm(n) * c(sqrt(var(mu_z2 + u2) / prop_zvar * (1-prop_zvar)))
z2 <- mu_z2 + u2 + e_z2

#### recover mediated effect x --(mediated)-> z ####
xz_coefs <- data.frame(summary(lm(z2 ~ 1 + x2))$coefficients[-1,])
xz_coefs_direct <- data.frame(summary(lm(z2 ~ 1 + y2 + x2))$coefficients[-c(1:2),])
xy_coefs <- data.frame(summary(lm(y2 ~ 1 + x2))$coefficients[-1,])
xz_coefs_mediated <- data.frame("Estimate" = xz_coefs$Estimate - xz_coefs_direct$Estimate,
                           "Std. Error" = sqrt(xz_coefs$Std..Error^2 + xz_coefs_direct$Std..Error^2)) #need to account for covariance
# hist(1 - abs(pt(xz_coefs_mediated$Estimate / xz_coefs_mediated$Std..Error, df = n-p) - 0.5) * 2)
# plot(xz_coefs$Estimate, xy_coefs$Estimate)
# plot(xz_coefs_mediated$Estimate, xy_coefs$Estimate)
summary(lm(xz_coefs_mediated$Estimate ~ xy_coefs$Estimate))$coefficients[-1,]

# plot(xz_coefs_mediated$Estimate, b_xy * b_xz)

#### SEM ####
colnames(x2) <- predictors <- paste0("x_", 1:p)

# Start constructing the SEM model
model <- ""
latent_effects <- "U =~ 1*y + 1*z"
direct_effects_to_y <- paste("y ~ cy * U + ", paste(paste0("b_xy_", 1:p, " * ", predictors), collapse = " + "))
direct_effects_to_z <- paste("z ~ cz * U + b_yz * y + ", paste(paste0("b_xz_", 1:p, " * ", predictors), 
                                          collapse = " + "))
model <- paste(latent_effects, direct_effects_to_y, direct_effects_to_z, sep = "\n")

indirect_effects <- paste0("indirect_", predictors, " := b_xy_", 1:p, " * b_yz")
total_effects <- paste0("total_", predictors, " := b_xz_", 1:p, " + (b_xy_", 1:p, " * b_yz)")
model <- paste(model, paste(indirect_effects, collapse = "\n"), 
               paste(total_effects, collapse = "\n"), sep = "\n")

# Create a data frame from the matrices/vectors
data <- data.frame(x2, y = y2, z = z2)

# Fit the SEM model
SEM_fit <- summary(sem(model, data = data))
SEM_fit$pe[SEM_fit$pe$label == "b_yz",]
summary(lm(xz_coefs$Estimate ~ xy_coefs$Estimate))$coefficients[-1,]

#### StanEM ####

# Prepare data for Stan
dat <- list(
  N = nrow(x2),
  P = ncol(x2),
  x = as.matrix(x2), # Convert x2 to matrix
  y = as.vector(y2), # Convert y2 to vector
  z = as.vector(z2)  # Convert z2 to vector
)

# Define Stan model as a string
stan_model <- "
data {
  int<lower=0> N; // number of observations
  int<lower=0> P; // number of predictors in x
  matrix[N, P] x; // matrix of predictors (instruments)
  vector[N] y;    // mediator
  vector[N] z;    // outcome
}

parameters {
  vector[P] b_x_y;    // coefficients for predictors on y
  real b_y_z;         // coefficient for y on z
  vector[P] b_x_z;    // coefficients for predictors on z
  real<lower=-1, upper=1> rho; // extra covariance between y and z not due to y -> z
  real<lower=0> sigma_y; // standard deviation for y
  real<lower=0> sigma_z; // standard deviation for z
}

model {
  // Priors
  b_x_y ~ std_normal(); // prior for coefficients of x on y
  b_y_z ~ std_normal(); // prior for coefficient of y on z
  b_x_z ~ std_normal(); // prior for coefficients of x on z
  rho ~ std_normal();   // prior for extra covariance (note: normal, truncated to -1 to 1)
  sigma_y ~ std_normal(); // prior for standard deviation of y
  sigma_z ~ std_normal(); // prior for standard deviation of z

  // Joint likelihood of y and z as bivariate normal
  for (n in 1:N) {
    vector[2] mu;
    cov_matrix[2] Sigma;
    
    // Means for y and z
    mu[1] = dot_product(b_x_y, x[n, ]); // mean of y given x
    mu[2] = b_y_z * y[n] + dot_product(b_x_z, x[n, ]); // mean of z given y and x

    // Covariance matrix for y and z
    Sigma[1, 1] = square(sigma_y);
    Sigma[2, 2] = square(sigma_z);
    Sigma[1, 2] = rho * sigma_y * sigma_z; // extra covariance term
    Sigma[2, 1] = Sigma[1, 2]; // symmetric covariance matrix

    // Bivariate normal distribution for y and z
    [y[n], z[n]] ~ multi_normal(mu, Sigma);
  }
}
"


mod <- cmdstan_model(write_stan_file(stan_model))
fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                      adapt_delta = 0.95, parallel_chains = 4,
                      refresh = 100, max_treedepth = 10, 
                      thin = 1, init = 0.1, data = dat)

summ <- fit$summary()
print(summ[order(summ$rhat, decreasing = T),])
samps <- data.frame(as_draws_df(fit$draws()))