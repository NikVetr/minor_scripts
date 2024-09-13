
data{
  
    int<lower = 1> p; // number of parameters per set
    int<lower = 1> n; // number of obs in original experiment
    matrix[p, 2] x_err;
    matrix<lower=0>[p, 2] sd_x_err;
    int<lower = 0> p0; // expected number of slab predictors in each condition
    
}

transformed data{
  
  // horseshoe hyperhyperparameters
  real<lower = 0> rho = sqrt(abs(1/mean(x_err)/(1-mean(x_err))));
  real<lower=1> pr = p * 1.0;
  real<lower = 0> p0r = p0 * 1.0;
  real<lower = 0> nr = n * 1.0;
  real<lower = 0> scale_global = p0r ./ (pr - p0r) / sqrt(nr); // scale for the half -t prior for tau
  real<lower = 1> nu_global = 1; // degrees of freedom for the half -t prior for tau
  real<lower = 1> nu_local = 1; // degrees of freedom for the half - t priors for lambdas
  real<lower = 0> slab_scale = 2; // slab scale for the regularized horseshoe
  real<lower = 0> slab_df = 4; // slab degrees of freedom for the regularized horseshoe
  
}

parameters{
  
  // MVN parameters
  array[p] vector[2] z;
  cholesky_factor_corr[2] L;
  
  // hierarchical horseshoe hyperparameters
  real<lower = 0> aux1_global;
  real<lower = 0> aux2_global;
  vector<lower = 0>[p] aux1_local;
  vector<lower = 0>[p] aux2_local;
  real<lower = 0> caux;
  
}

transformed parameters{
  
  // horseshoe hyperparameter transformations
  real<lower = 0> tau = aux1_global .* sqrt(aux2_global) .* scale_global; // global shrinkage parameter
  real<lower = 0> c = slab_scale .* sqrt(caux); // slab scale
  
  vector<lower = 0>[p] lambda =  aux1_local .* sqrt(aux2_local);; // local shrinkage parameter
  vector<lower = 0>[p] lambda_tilde =  sqrt(c^2 * lambda^2 ./ (c^2 + tau^2 * lambda^2));; // truncated local shrinkage parameter
  array[p] vector[2] x; // latent effect sizes
  
  //combining everything together
  for(i in 1:2){
    x[,i] =  to_array_1d(to_vector(z[,i]) .* lambda_tilde * tau);
  }
  
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
  L ~ lkj_corr_cholesky(1);
  z ~ multi_normal_cholesky(rep_vector(0, 2), L);

  // observed effect model
  // x_err[,1] ./ sd_x_err[,1] - to_vector(x[,1]) ~ std_normal();
  // x_err[,2] ./ sd_x_err[,2] - to_vector(x[,2]) ~ std_normal();
  x_err[,1] ~ normal(x[,1], sd_x_err[,1]);
  x_err[,2] ~ normal(x[,2], sd_x_err[,2]);
  
}

generated quantities{
  
  // recover correlation coefficient
  real r = multiply_lower_tri_self_transpose(L)[1,2];
  
}
