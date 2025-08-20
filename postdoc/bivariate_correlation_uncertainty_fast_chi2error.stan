data{
  int<lower=1>           p;
  matrix[p,2]            x_err;
  matrix<lower=0>[p,2]   sd_x_err;           // reported s.e.’s
  array[p,2] int<lower=1> df;                // same df as before
}
parameters{
  // latent coefficients (non-centred)
  vector[p]              z1, z2;
  vector<lower=0>[2]     sigma_x;
  real<lower=-1,upper=1> r;

  // put σ_e on **log** scale for better geometry
  vector[p]              log_sigma_e1;
  vector[p]              log_sigma_e2;
}
transformed parameters{
  vector<lower=0>[p] sigma_e1 = exp(log_sigma_e1);
  vector<lower=0>[p] sigma_e2 = exp(log_sigma_e2);

  matrix[p,2] x;
  x[,1] = sigma_x[1] .* z1;
  x[,2] = sigma_x[2] .* ( r * z1 + sqrt(1 - square(r)) .* z2 );
}
model{
  /* priors */
  sigma_x      ~ normal(0, 2);
  z1           ~ std_normal();
  z2           ~ std_normal();

  // half-Normal centred on the reported s.e.
  log_sigma_e1 ~ normal(log(sd_x_err[,1]), 5);   // sd ≈ factor ≈ e
  log_sigma_e2 ~ normal(log(sd_x_err[,2]), 5);

  // pt estimates
  x_err[,1] ~ normal(x[,1], sigma_e1);
  x_err[,2] ~ normal(x[,2], sigma_e2);

  // SEs
  vector[p] nu1 = to_vector(df[,1]);
  vector[p] nu2 = to_vector(df[,2]);
  vector[p] s2_1 = square(sd_x_err[,1]);
  vector[p] s2_2 = square(sd_x_err[,2]);
  
  target += chi_square_lpdf( nu1 .* s2_1 ./ square(sigma_e1) | nu1 )
            + sum( log(nu1) - 2 * log_sigma_e1 );
  target += chi_square_lpdf( nu2 .* s2_2 ./ square(sigma_e2) | nu2 )
            + sum( log(nu2) - 2 * log_sigma_e2 );

}


