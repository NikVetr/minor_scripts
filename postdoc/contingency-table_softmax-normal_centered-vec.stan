data {
  int<lower=1> k; // number of rows
  int<lower=1> p; // number of columns
  array[k*p] int<lower=1, upper=k> row_id; // row indices
  array[k*p] int<lower=1, upper=p> col_id; // column indices
  array[k*p] int<lower=1, upper=k*p> cell_id; // cell indices
  array[k*p] int<lower=0> count;
  real<lower=0> gamma_sd_sd; // standard deviation factor for the interactions scale
}

transformed data {
  int<lower=1> kp = k * p;
}

parameters {
  // Use sum_to_zero_vector for the base effects
  sum_to_zero_vector[k] alpha_base; // base row effects (sum to zero)
  real<lower=0> alpha_sd;          // scale for row effects
  sum_to_zero_vector[p] beta_base;  // base col effects (sum to zero)
  real<lower=0> beta_sd;           // scale for col effects
  sum_to_zero_vector[kp] gamma_base; // base interaction effects (sum to zero)
  real<lower=0> gamma_sd;          // scale for interaction effects
}

transformed parameters {
  // Construct the scaled effects; these will also sum to zero
  vector[k] alpha_centered = alpha_base * alpha_sd; // row effects
  vector[p] beta_centered = beta_base * beta_sd;   // column effects
  vector[kp] gamma_centered = gamma_base * gamma_sd * gamma_sd_sd; // interaction effects

  // Construct latents using the centered effects
  vector[kp] eta = alpha_centered[row_id] + beta_centered[col_id] + gamma_centered[cell_id];
}

model {
  // Priors
  alpha_base ~ std_normal();
  beta_base ~ std_normal();
  gamma_base ~ std_normal();

  alpha_sd ~ std_normal();
  beta_sd ~ std_normal();
  gamma_sd ~ std_normal();

  // Likelihood
  count ~ multinomial_logit(eta);
}

