
data{
  
    int<lower = 1> p; // number of parameters per set
    matrix[p, 2] x_err;
    matrix<lower=0>[p, 2] sd_x_err;
    array[2] int<lower = 0> p0; // expected number of slab predictors in each condition
    
}

transformed data{
  
  // horseshoe hyperhyperparameters
  vector<lower = 0>[2] rho;
  rho[1] = sqrt(abs(1/mean(x_err[,1])/(1-mean(x_err[,1]))));
  rho[2] = sqrt(abs(1/mean(x_err[,2])/(1-mean(x_err[,2]))));
  
  real<lower=1> pr = p * 1.0;
  vector<lower = 0>[2] p0v = to_vector(p0);
  vector<lower = 0>[2] scale_global = p0v ./ (pr - p0v) .* rho / sqrt(pr); // scale for the half -t prior for tau
  vector<lower = 1>[2] nu_global = rep_vector(1, 2); // degrees of freedom for the half -t prior for tau
  vector<lower = 1>[2] nu_local = rep_vector(1, 2); // degrees of freedom for the half - t priors for lambdas
  vector<lower = 0>[2] slab_scale = rep_vector(2, 2); // slab scale for the regularized horseshoe
  vector<lower = 0>[2] slab_df = rep_vector(4, 2); // slab degrees of freedom for the regularized horseshoe
  
}

parameters{
  
  // MVN parameters
  array[p] vector[2] z;
  cholesky_factor_corr[2] L;
  
  // horseshoe hyperparameters
  vector<lower = 0>[2] aux1_global;
  vector<lower = 0>[2] aux2_global;
  array[2] vector<lower = 0>[p] aux1_local;
  array[2] vector<lower = 0>[p] aux2_local;
  vector<lower = 0>[2] caux;
  
}

transformed parameters{
  
  // horseshoe hyperparameter transformations
  vector<lower = 0>[2] tau = aux1_global .* sqrt(aux2_global) .* scale_global; // global shrinkage parameter
  vector<lower = 0>[2] c = slab_scale .* sqrt(caux); // slab scale
  
  array[2] vector<lower = 0>[p] lambda; // local shrinkage parameter
  array[2] vector<lower = 0>[p] lambda_tilde; // truncated local shrinkage parameter
  array[p] vector[2] x; // latent effect sizes
  for(i in 1:2){
    lambda[i] =  aux1_local[i] .* sqrt(aux2_local[i]);
    lambda_tilde[i] =  sqrt(c[i]^2 * lambda[i]^2 ./ (c[i]^2 + tau[i]^2 * lambda[i]^2));
    x[,i] =  to_array_1d(to_vector(z[,i]) .* lambda_tilde[i] * tau[i]);
  }
  
}

model{
  
  // latent effect model (horseshoe)
  // half -t priors for lambdas and tau, and inverse - gamma for c^2
  for(i in 1:2){
    aux1_local[i] ~ normal(0, 1);
    aux2_local[i] ~ inv_gamma(0.5 * nu_local[i], 0.5 * nu_local[i]);
    aux1_global[i] ~ normal(0, 1);
    aux2_global[i] ~ inv_gamma(0.5 * nu_global[i], 0.5 * nu_global[i]);  
    caux[i] ~ inv_gamma(0.5 * slab_df[i], 0.5 * slab_df[i]);
  }
  
  // observed effect model
  // x_err[,1] ./ sd_x_err[,1] - to_vector(x[,1]) ~ std_normal();
  // x_err[,2] ./ sd_x_err[,2] - to_vector(x[,2]) ~ std_normal();
  L ~ lkj_corr_cholesky(1);
  z ~ multi_normal_cholesky(rep_vector(0, 2), L);
  
  x_err[,1] ~ normal(x[,1], sd_x_err[,1]);
  x_err[,2] ~ normal(x[,2], sd_x_err[,2]);

}

generated quantities{
  
  // recover correlation coefficient
  real r = multiply_lower_tri_self_transpose(L)[1,2];
  
}
