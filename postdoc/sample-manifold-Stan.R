library(cmdstanr)
library(posterior)
library(ggplot2)

# Define the Stan model as a string
stan_model_code <- "
data {
  int<lower=1> N;  // Number of samples
}
parameters {
  vector[4] p;  // (x, y, z, w) coordinates in 4D
}
model {
  // Constraint 1: Must be on the unit sphere
  target += -1000 * abs(dot_self(p) - 1);

  // Constraint 2: Must satisfy the hyperellipsoid equation
  target += -1000 * abs(2 * p[1]^2 + p[2]^2 + 0.75 * p[3]^2 + 0.25 * p[4]^2 - 1);

  // Constraint 3: Must satisfy the hyperplane equation
  target += -1000 * abs(0.52 * p[1] - 0.59 * p[2] + 0.15 * p[3] - 0.60 * p[4]);

  // Implicit uniform prior (HMC ensures uniform exploration)
}
"

# Write the Stan model to a temporary file
stan_file <- tempfile(fileext = ".stan")
writeLines(stan_model_code, stan_file)

# Compile the Stan model
mod <- cmdstan_model(stan_file)

# Define the number of samples
N_samples <- 1000

# Run HMC sampling
fit <- mod$sample(
  data = list(N = N_samples),
  iter_sampling = 2000,  # Number of post-warmup samples
  iter_warmup = 1000,     # Warmup iterations
  chains = 4,             # Number of chains
  parallel_chains = 4     # Run in parallel
)

# Extract samples as a data frame
samples <- fit$draws(format = "df")

# Extract (x, y, z, w) columns
samps <- data.frame(samples[, c("p[1]", "p[2]", "p[3]", "p[4]")])
samps <- cbind(samps, samps^2)
colnames(samps) <- c(c("x", "y", "z", "w"), paste0(c("x", "y", "z", "w"), 2))
samps <- cbind(samps, samps^2)
colnames(samps)[5:8] <- colnames(samps)[5:8]

# Print a summary of the sampled points
summary(lm(y ~ 1 + z2 + w2 + z * w + z + w, samps))
print(summary(samps))

# Basic visualization of (z, w) projections
ggplot(samps, aes(x = z, y = w)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Projection of Samples onto (z, w) Plane",
       x = "z", y = "w")
