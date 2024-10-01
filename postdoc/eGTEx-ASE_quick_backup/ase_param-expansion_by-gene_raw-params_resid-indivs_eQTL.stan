data {
  int<lower=0> n;
  int<lower=0> n_indiv;
  array[n] int<lower=1, upper=n_indiv> indiv_j;
  array[n_indiv] int<lower=1, upper=2> eQTL_het;
  array[n] int<lower=0> count;
  array[n] int<lower=0> total;
  vector<lower=0.5, upper=1>[2] bounds;
}

transformed data {
  real<lower=0, upper=1> slab_width = bounds[2] - bounds[1];
  real<lower=0, upper=1> AB_width = bounds[1] - 0.5;
  real<lower=0, upper=1> loh_width = 1 - bounds[2];
}

parameters {
  
  // slab parameters //
  // for truncated normal //
  
  //symmetric location parameter
  real<lower=bounds[1], upper=bounds[2]> mu_base;
  real<lower=-slab_width, upper=slab_width> mu_eQTL_diff;
  
  //sigma parameter
  real<lower=-10, upper=10> sigma_base;
  real<lower=-10, upper=10> sigma_eQTL_diff;
  
  vector<lower=bounds[1], upper=bounds[2]>[n_indiv] theta_slab;

  //logit mixture param for allelic balance
  real logodds_AB_base;
  real logodds_AB_eQTL_diff;
  
  //laplace probability params for AB 
  vector<lower=1-bounds[1], upper=bounds[1]>[n_indiv] theta_AB;
  real<lower=0, upper=2> sigma_theta_AB;
  
  //logit mixture param for loh (vs slab)
  real cond_logodds_loh_base;
  real cond_logodds_loh_eQTL_diff;
  
  //exponential probability params for loh 
  vector<lower=0, upper=1-bounds[2]>[n_indiv] theta_loh;
  real<lower=0, upper=1> sigma_theta_loh;
  
}

transformed parameters {
  
  //////////////////////////
  // slab location params //
  //////////////////////////
  
  //construct positive location vector
  vector[2] mu = mu_base + to_vector([-mu_eQTL_diff/2, mu_eQTL_diff/2]);
  
  ////////////////////////
  // slab scale params //
  ////////////////////////
  
  //construct positive sigmailon scale vector
  vector[2] sigma = exp(sigma_base + to_vector([-sigma_eQTL_diff/2, sigma_eQTL_diff/2]));
  
  ////////////////////
  // MIXTURE PARAMS //
  ////////////////////
  
  //get mixing weightability
  vector[2] prob_AB_vec = inv_logit(logodds_AB_base + 
    to_vector([-logodds_AB_eQTL_diff/2, logodds_AB_eQTL_diff/2]));
    
  vector[2] cond_prob_loh_vec = inv_logit(cond_logodds_loh_base +
    to_vector([-cond_logodds_loh_eQTL_diff/2, cond_logodds_loh_eQTL_diff/2]));
  
  //find each of the mixture probabilities
  array[2, 5] real mix_prob;
  for (i in 1:2) {
    mix_prob[i,1] = prob_AB_vec[i]; //probability in allelic balance
    mix_prob[i,2] = ((1-prob_AB_vec[i]) * cond_prob_loh_vec[i]) / 2; //probability not in AB but in LoH (tail 1)
    mix_prob[i,3] = mix_prob[i,2]; //probability not in AB but in LoH (tail 2)
    mix_prob[i,4] = ((1-prob_AB_vec[i]) * (1-cond_prob_loh_vec[i])) / 2; //probability not in AB but in slab (tail 1)
    mix_prob[i,5] = mix_prob[i,4]; //probability not in AB but in slab (tail 2)
  }
  array[2, 5] real log_mix_prob = log(mix_prob);
  
}

model {
  
  //////////////////////////// 
  // priors and hyperpriors //
  //////////////////////////// 
  
  mu_base ~ normal((bounds[1] + bounds[2])/2, slab_width/6);
  mu_eQTL_diff ~ normal(0, slab_width/4);
  
  sigma_base ~ std_normal();
  sigma_eQTL_diff ~ std_normal();
  
  logodds_AB_base ~ std_normal();
  logodds_AB_eQTL_diff ~ normal(0, 0.5);
  
  cond_logodds_loh_base ~ std_normal();
  cond_logodds_loh_eQTL_diff ~ normal(0, 0.5);
  
  sigma_theta_loh ~ normal(0, 0.5);
  sigma_theta_AB ~ normal(0, 0.5);
  
  theta_AB ~ double_exponential(0.5, sigma_theta_AB);
  theta_slab ~ normal(mu[eQTL_het], sigma[eQTL_het]);
  theta_loh ~ exponential(1/sigma_theta_loh);
  
  ////////////////////////
  // Mixture Likelihood //
  ////////////////////////
  
  //average binomial probabilities across mixture
  for (i in 1:n) {
    array[5] real lp;
    int i_het = eQTL_het[indiv_j[i]];
    lp[1] = log_mix_prob[i_het,1] + binomial_lpmf(count[i] | total[i], theta_AB[indiv_j[i]]);
    lp[2] = log_mix_prob[i_het,2] + binomial_lpmf(count[i] | total[i], theta_loh[indiv_j[i]]);
    lp[3] = log_mix_prob[i_het,3] + binomial_lpmf(count[i] | total[i], 1-theta_loh[indiv_j[i]]);
    lp[4] = log_mix_prob[i_het,4] + binomial_lpmf(count[i] | total[i], theta_slab[indiv_j[i]]);
    lp[5] = log_mix_prob[i_het,5] + binomial_lpmf(count[i] | total[i], 1-theta_slab[indiv_j[i]]);
    target += log_sum_exp(lp);
  }
  
}
