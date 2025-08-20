data {
  int<lower = 0> n; // number of observations
  int<lower = 0> p; // number of predictors
  vector[n] y; // outputs
  matrix[n, p] x; // inputs
  real<lower = 0> scale_icept; // prior std for the intercept
  real<lower = 0> scale_global; // scale for the half -t prior for tau
  real<lower = 1> nu_global; // degrees of freedom for the half -t prior for tau
  real<lower = 1> nu_local; // degrees of freedom for the half - t priors for lambdas
  real<lower = 0> slab_scale; // slab scale for the regularized horseshoe
  real<lower = 0> slab_df; // slab degrees of freedom for the regularized horseshoe
}

parameters {
  real logsigma;
  real beta0;
  vector[p] z;
  real<lower = 0> aux1_global;
  real<lower = 0> aux2_global;
  vector<lower = 0>[p] aux1_local;
  vector<lower = 0>[p] aux2_local;
  real<lower = 0> caux;
}

transformed parameters {
  real<lower = 0> sigma; // noise std
  real<lower = 0> tau; // global shrinkage parameter
  vector<lower = 0>[p] lambda; // local shrinkage parameter
  vector<lower = 0>[p] lambda_tilde; // ’ truncated ’ local shrinkage parameter
  real<lower = 0> c; // slab scale
  vector[p] beta; // regression coefficients
  vector[n] f; // latent function values
  sigma =  exp(logsigma);
  lambda =  aux1_local .* sqrt(aux2_local);
  tau =  aux1_global * sqrt(aux2_global) * scale_global * sigma;
  c =  slab_scale * sqrt(caux);
  lambda_tilde =  sqrt(c^2 * lambda^2 ./(c^2 + tau^2* lambda^2));
  beta =  z .* lambda_tilde * tau;
  f =  beta0 + x * beta;
}

model {
  // half -t priors for lambdas and tau, and inverse - gamma for c^2
  logsigma ~ normal(0, 1000);
  z ~ normal(0, 1);
  aux1_local ~ normal(0, 1);
  aux2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local);
  aux1_global ~ normal(0, 1);
  aux2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  beta0 ~ normal(0, scale_icept);
  y ~ normal(f, sigma);
}
