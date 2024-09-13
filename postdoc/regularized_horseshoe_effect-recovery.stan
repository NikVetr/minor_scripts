
data{
  
    int<lower = 1> p; // number of parameters per set
    int<lower = 1> n; // number of obs in the original inference
    vector[p] x_err;
    vector<lower=0>[p] sd_x_err;
    int<lower = 0> p0; // expected number of slab predictors in each condition
    
}

transformed data{
  
  // horseshoe hyperhyperparameters
  real<lower = 0> pseudo_sigma = sqrt(abs(1/mean(x_err)/(1-mean(x_err))));

  real<lower=1> pr = p * 1.0;
  real<lower = 0> p0r = p0 * 1.0;
  real<lower = 0> scale_global = p0r / (pr - p0r) / sqrt(n); // scale for the half -t prior for tau
  real<lower = 1> nu_global = 1; // degrees of freedom for the half -t prior for tau
  real<lower = 1> nu_local = 1; // degrees of freedom for the half - t priors for lambdas
  real<lower = 0> slab_scale = 2; // slab scale for the regularized horseshoe
  real<lower = 0> slab_df = 4; // slab degrees of freedom for the regularized horseshoe
  
}

parameters{
  
  // MVN parameters
  vector[p] z;
  
  // hierarchical horseshoe hyperparameters
  real logsigma;
  real<lower = 0> aux1_global;
  real<lower = 0> aux2_global;
  vector<lower = 0>[p] aux1_local;
  vector<lower = 0>[p] aux2_local;
  real<lower = 0> caux;
  
}

transformed parameters{
  
  // horseshoe hyperparameter transformations
  real <lower =0> sigma = exp(logsigma);
  real<lower = 0> tau = aux1_global * sqrt(aux2_global) * scale_global * sigma; // global shrinkage parameter
  real<lower = 0> c = slab_scale * sqrt(caux); // slab scale
  vector<lower = 0>[p] lambda =  aux1_local .* sqrt(aux2_local); // local shrinkage parameter
  vector<lower = 0>[p] lambda_tilde=  sqrt(c^2 * square(lambda) ./ (c^2 + tau^2 * square(lambda))); // truncated local shrinkage parameter
  vector[p] x =  z .* lambda_tilde * tau; // latent effect sizes
  
}

model{
  
  // latent effect model (horseshoe)
  // half -t priors for lambdas and tau, and inverse - gamma for c^2
  aux1_local ~ std_normal();
  aux2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local);
  aux1_global ~ std_normal();
  aux2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global);  
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  // priors on latent effects
  z ~ std_normal();

  // observed effect model
  x_err ~ normal(x, sd_x_err);
  
}
