data{
  
    int<lower = 1> p; // number of parameters per set
    matrix[p, 2] x_err;
    matrix<lower=0>[p, 2] sd_x_err;
    int<lower = 0> p0; // expected number of non-null parameters
    
}

parameters{
  
  // MVN parameters
  array[p] vector[2] z;
  cholesky_factor_corr[2] L;
  vector<lower=0>[2] sigma_x;
  
  // hierarchical sparsity parameters
  real mu_resc;
  vector[p] logit_resc;
  
}

transformed parameters{
  
  //construct rescale parameter in (0,1)
  // vector[p] resc = 1 / (1 + exp(-100 * (inv_logit(mu_resc + logit_resc) - 0.5)));
  vector[p] resc = inv_logit((mu_resc + logit_resc)*100);
  array[p] vector[2] x; // latent effect sizes
  
  //combining everything together
  for(i in 1:2){
    x[,i] =  to_array_1d(to_vector(z[,i]) .* resc * sigma_x[i]);
  }
  
}

model{
  
  // priors on sparsity params
  mu_resc ~ std_normal();
  logit_resc ~ std_normal();
  
  // priors on latent effects
  L ~ lkj_corr_cholesky(1);
  z ~ multi_normal_cholesky(rep_vector(0, 2), L);
  sigma_x ~ std_normal();

  // observed effect model
  x_err[,1] ~ normal(x[,1], sd_x_err[,1]);
  x_err[,2] ~ normal(x[,2], sd_x_err[,2]);
  
}

generated quantities{
  
  // recover correlation coefficient
  real r = multiply_lower_tri_self_transpose(L)[1,2];
  
}
