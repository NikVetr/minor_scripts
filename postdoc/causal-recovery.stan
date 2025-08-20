data  {
  int<lower=1>  p;
  vector[p]     b_xy_est, b_xy_se;
  vector[p]     b_xz_est, b_xz_se;
}
parameters {
  vector[p]     b_xy_raw;
  real          beta_yz;
  real          alpha;
  real<lower=0> sigma_xy;
  real<lower=0> tau;           // √τ  = pleiotropy SD
}
transformed parameters {
  vector[p] b_xy = b_xy_raw * sigma_xy;
  vector[p] mu_xz = alpha + beta_yz * b_xy;        // mean relation
  vector[p] var_xz = square(b_xz_se) + square(tau);// SE²  + τ²
}
model {
  // priors
  beta_yz    ~ normal(0,5);
  alpha      ~ normal(0,5);
  b_xy_raw   ~ std_normal();
  sigma_xy   ~ student_t(3,0,1);   // weakly-informative half-t
  tau        ~ student_t(3,0,1);

  // first-stage error propagation
  b_xy_est   ~ normal(b_xy, b_xy_se);

  // second stage with marginalised δ_k
  b_xz_est   ~ normal(mu_xz, sqrt(var_xz));
}
