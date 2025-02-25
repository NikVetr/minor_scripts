# Load necessary library
library(cmdstanr)

# Step 1: Define a toy Stan model and compile it
stan_code <- "
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real theta;
}
model {
  mu ~ normal(0, 1);
  theta ~ normal(0, 1);
  y ~ normal(mu + theta, 1);
}
"
# Save and compile the model
file <- write_stan_file(stan_code)
mod <- cmdstan_model(file)

# Step 2: Generate some toy data
N <- 100
y <- rnorm(N, mean = 0, sd = 1)
data <- list(N = N, y = y)

# Step 3: Run 4 independent chains
fit_1 <- mod$sample(data = data, chains = 1, seed = 1)
fit_2 <- mod$sample(data = data, chains = 1, seed = 2)
fit_3 <- mod$sample(data = data, chains = 1, seed = 3)
fit_4 <- mod$sample(data = data, chains = 1, seed = 4)

# Step 4: Combine the CSV outputs from all chains into a single CmdStanMCMC object
output_files <- c(fit_1$output_files(), fit_2$output_files(), fit_3$output_files(), fit_4$output_files())
combined_csv <- cmdstanr::read_cmdstan_csv(output_files)

# Step 5: Extract draws and summarize them
draws_array <- combined_csv$post_warmup_draws
draws <- posterior::as_draws_array(draws_array)
summ <- posterior::summarise_draws(draws)

#compare to individual chains
summ
fit_1$summary()
fit_2$summary()
fit_3$summary()
fit_4$summary()
