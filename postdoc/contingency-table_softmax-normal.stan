functions {
  array[] int drop_index(array[] int x, int k) {
    int n = size(x);
    array[n - 1] int out;
    int j = 1;
    for (i in 1:n) {
      if (i != k) {
        out[j] = x[i];
        j += 1;
      }
    }
    return out;
  }
}

data {
  int<lower=1> k; // number of rows
  int<lower=1> p; // number of columns
  array[k*p] int<lower=1, upper=k> row_id; // row indices 
  array[k*p] int<lower=1, upper=p> col_id; // column indices
  array[k*p] int<lower=1, upper=k*p> cell_id; // cell indices
  array[k*p] int<lower=0> count;
  int<lower=1, upper=k> ref_row; // reference for 0-row effect
  int<lower=1, upper=p> ref_col; // reference for 0-column effect
  int<lower=1, upper=k*p> ref_cell; // reference for 0-interaction effect
  real<lower=0> gamma_sd_sd; // standard deviation for the interactions cale
}

transformed data {
  int<lower=1> kp = k * p;

  array[k] int row_inds;
  for (i in 1:k) row_inds[i] = i;

  array[p] int col_inds;
  for (i in 1:p) col_inds[i] = i;

  array[kp] int cell_inds;
  for (i in 1:kp) cell_inds[i] = i;

  array[k-1] int nref_row = drop_index(row_inds, ref_row);
  array[p-1] int nref_col = drop_index(col_inds, ref_col);
  array[kp-1] int nref_cell = drop_index(cell_inds, ref_cell);
}

parameters {
  vector[k - 1] alpha_raw; // row intercepts
  real<lower=0> alpha_sd;
  vector[p - 1] beta_raw; // col intercepts
  real<lower=0> beta_sd;
  vector[kp - 1] gamma_raw; // interaction effects
  real<lower=0> gamma_sd;
}

transformed parameters {
  vector[k] alpha; // row effects
  alpha[ref_row] = 0; // fix "largest" entry for identifiability
  alpha[nref_row] = alpha_raw * alpha_sd;
  
  vector[p] beta; // column effects
  beta[ref_col] = 0; // fix "largest" entry for identifiability
  beta[nref_col] = beta_raw * beta_sd;
  
  vector[kp] gamma; // interaction effects
  gamma[ref_cell] = 0; // fix "largest" entry for identifiability
  gamma[nref_cell] = gamma_raw * gamma_sd * gamma_sd_sd;
  
  vector[kp] eta = alpha[row_id] + beta[col_id] + gamma[cell_id]; // construct latents
}

model {
  // priors
  alpha_raw ~ std_normal();
  beta_raw ~ std_normal();
  gamma_raw ~ std_normal();
  alpha_sd ~ std_normal();
  beta_sd ~ std_normal();
  gamma_sd ~ std_normal();
  
  // likelihood
  count ~ multinomial_logit(eta);
}

generated quantities {
  // recover centered parameters (to accommodate fixing for identifiability)
  vector[k] alpha_centered = alpha - mean(alpha);
  vector[p] beta_centered = beta - mean(beta);
  vector[kp] gamma_centered = gamma - mean(gamma);
}
