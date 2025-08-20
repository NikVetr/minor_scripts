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
  vector[k] alpha;
  vector[p] beta;
  array[k*p] int<lower=1, upper=k> row_id; // row indices 
  array[k*p] int<lower=1, upper=p> col_id; // column indices
  array[k*p] int<lower=1, upper=k*p> cell_id; // cell indices
  array[k*p] int<lower=0> count;
  int<lower=1, upper=k*p> ref_cell; // reference for 0-interaction effect
}

transformed data {
  int<lower=1> kp = k * p;
  array[kp] int cell_inds;
  for (i in 1:kp) cell_inds[i] = i;
  array[kp-1] int nref_cell = drop_index(cell_inds, ref_cell);
}

parameters {
  vector[kp - 1] gamma_raw; // interaction effects
  real<lower=0> gamma_sd;
}

transformed parameters {
  vector[kp] gamma; // interaction effects
  gamma[ref_cell] = 0; // fix "largest" entry for identifiability
  gamma[nref_cell] = gamma_raw * gamma_sd;
  vector[kp] eta = alpha[row_id] + beta[col_id] + gamma[cell_id]; // construct latents
}

model {
  // priors
  gamma_raw ~ std_normal();
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
