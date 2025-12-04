data{
  int<lower=1> n;
  int<lower=1> d;
  matrix[n, d] x;
  vector[n] y;                     // ignored when prior_only = 1
  int<lower=1> G;
  array[d] int<lower=1, upper=G> group;

  vector<lower=0>[G] hs_tau0;
  vector<lower=0>[G] slab_scale;
  real<lower=1> slab_df;

  int<lower=0, upper=1> prior_only;  // NEW
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
  lambda ~ cauchy(0, 1);
  tau ~ student_t(3, 0, hs_tau0);

  for (g in 1:G)
    c2[g] ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df * square(slab_scale[g]));

  a     ~ normal(0, 1);
  sigma ~ student_t(3, 0, 2);
  df    ~ student_t(3, 2, 5);

  // likelihood guarded by switch
  if (prior_only == 0)
    y ~ student_t(df, a + x * beta, sigma);
}
generated quantities{
  // shrinkage diagnostics
  array[d] real kappa;
  real p_eff = 0;
  for (j in 1:d) {
    int g = group[j];
    real num = tau[g]*lambda[j];
    real den = sqrt( c2[g] / lam2[j] ) + num;
    kappa[j] = 1 - square(num/den);
    p_eff += 1 - kappa[j];
  }

  // prior predictive draws for y
  vector[n] mu = a + x * beta;
  vector[n] y_prior;
  for (i in 1:n)
    y_prior[i] = student_t_rng(df, mu[i], sigma);
}
