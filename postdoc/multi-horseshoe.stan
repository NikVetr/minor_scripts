data{
  int<lower=1> n;
  int<lower=1> d;
  matrix[n, d] x;
  vector[n] y;
  int<lower=1> G;
  array[d] int<lower=1, upper=G> group;

  vector<lower=0>[G] hs_tau0;
  vector<lower=0>[G] slab_scale;
  real<lower=1> slab_df;
}
parameters{
  real a;
  real<lower=0> sigma;
  real<lower=2> df;

  vector[d] z_beta;
  vector<lower=0>[d] lambda;
  vector<lower=0>[G] tau;
  vector<lower=0>[G] c2;
}
transformed parameters{
  vector[d] beta;
  vector[d] lam2 = square(lambda);
  vector[G] tau2 = square(tau);

  for (j in 1:d) {
    int g = group[j];
    real lam_tilde = sqrt( c2[g] * lam2[j] / (c2[g] + tau2[g] * lam2[j]) );
    beta[j] = z_beta[j] * tau[g] * lam_tilde;
  }
}
model{
  // priors
  z_beta ~ normal(0, 1);

  // half-Cauchy via <lower=0> and Cauchy(0, scale)
  lambda ~ cauchy(0, 1);
  tau ~ student_t(3, 0, hs_tau0);

  // slab: c2 ~ Inv-Gamma(ν/2, ν s^2/2)
  for (g in 1:G)
    c2[g] ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df * square(slab_scale[g]));

  a     ~ normal(0, 1);     // y is z-scored; tighter is fine
  sigma ~ student_t(3, 0, 2);
  df ~ student_t(3, 2, 5);
  
  y ~ student_t(df, a + x * beta, sigma);
}
generated quantities{
  // shrinkage diagnostics (effective number of nonzeros)
  array[d] real kappa;  // 1 = fully shrunk, 0 = unshrunk
  real p_eff = 0;
  for (j in 1:d) {
    int g = group[j];
    real num = tau[g]*lambda[j];
    real den = sqrt( c2[g] / lam2[j] ) + num;  // from RS horseshoe algebra
    kappa[j] = 1 - square(num/den);            // ≈ shrinkage factor
    p_eff += 1 - kappa[j];
  }
}
