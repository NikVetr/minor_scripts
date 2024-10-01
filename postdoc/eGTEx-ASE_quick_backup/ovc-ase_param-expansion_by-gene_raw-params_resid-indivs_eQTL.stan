data {
  int<lower=0> n;
  int<lower=0> n_indiv;
  int<lower=0> n_grade;
  array[n] int<lower=1, upper=n_indiv> indiv_j;
  array[n] int<lower=1, upper=2> ovc;
  array[n_indiv] int<lower=1, upper=n_grade> grade_k;
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
  real<lower=-slab_width, upper=slab_width> mu_ovc;
  real<lower=-slab_width, upper=slab_width> mean_mu_grade;
  vector<lower=-slab_width, upper=slab_width>[n_grade-2] mu_grade;
  real<lower=0, upper=0.5> sd_mu_grade;
  
  //sigma parameter
  real<lower=-10, upper=10> sigma_base;
  real<lower=-10, upper=10> sigma_ovc;
  real<lower=-10, upper=10> mean_sigma_grade;
  vector<lower=-10, upper=10>[n_grade-2] sigma_grade;
  real<lower=0, upper=0.5> sd_sigma_grade;
  
  vector<lower=bounds[1], upper=bounds[2]>[n_indiv] theta_slab;

  //logit mixture param for allelic balance
  real logodds_AB_base;
  real<lower=-10, upper=10> logodds_AB_ovc;
  real<lower=-10, upper=10> mean_logodds_AB_grade;
  vector<lower=-10, upper=10>[n_grade-2] logodds_AB_grade;
  real<lower=0, upper=0.5> sd_logodds_AB_grade;
  
  //laplace probability params for AB 
  vector<lower=1-bounds[1], upper=bounds[1]>[n_indiv] theta_AB;
  real<lower=0, upper=2> sigma_theta_AB;
  
  //logit mixture param for loh (vs slab)
  real cond_logodds_loh_base;
  real<lower=-10, upper=10> cond_logodds_loh_ovc;
  real<lower=-10, upper=10> mean_cond_logodds_loh_grade;
  vector<lower=-10, upper=10>[n_grade-2] cond_logodds_loh_grade;
  real<lower=0, upper=0.5> sd_cond_logodds_loh_grade;
  
  //exponential probability params for loh 
  vector<lower=0, upper=1-bounds[2]>[n_indiv] theta_loh;
  real<lower=0, upper=1> sigma_theta_loh;
  
}

transformed parameters {
  
  //////////////////////////
  // slab location params //
  //////////////////////////
  
  //construct cumulative grade effect
  vector[n_grade] mu_grade_csum;
  mu_grade_csum[1] = -mu_ovc/2;
  mu_grade_csum[2] = mu_ovc/2;
  for (i in 3:n_grade) {
    mu_grade_csum[i] = mu_grade_csum[i-1] + mu_grade[i-2] * sd_mu_grade + mean_mu_grade;
  }
  
  //construct positive location vector
  vector[n_grade] mu = mu_base + mu_grade_csum;
  
  ////////////////////////
  // slab scale params //
  ////////////////////////
  
  //construct grade effect
  vector[n_grade] sigma_grade_csum;
  sigma_grade_csum[1] = -sigma_ovc/2;
  sigma_grade_csum[2] = sigma_ovc/2;
  for (i in 3:n_grade) {
    sigma_grade_csum[i] = sigma_grade_csum[i-1] + sigma_grade[i-2] * sd_sigma_grade + mean_sigma_grade;
  }
  
  //construct positive scale vector
  vector[n_grade] sigma = exp(sigma_base + sigma_grade_csum);
  
  ////////////////////
  // MIXTURE PARAMS //
  ////////////////////
  
  //for being in the AB state
  vector[n_grade] logodds_AB_grade_csum;
  logodds_AB_grade_csum[1] = -logodds_AB_ovc/2;
  logodds_AB_grade_csum[2] = logodds_AB_ovc/2;
  for (i in 3:n_grade) {
    logodds_AB_grade_csum[i] = logodds_AB_grade_csum[i-1] + logodds_AB_grade[i-2] * sd_logodds_AB_grade + mean_logodds_AB_grade;
  }
  
  //construct positive probability vector
  vector[n_grade] prob_AB = inv_logit(logodds_AB_base + logodds_AB_grade_csum);
  
  
  //for being in the loh state when not in the AB state
  vector[n_grade] cond_logodds_loh_grade_csum;
  cond_logodds_loh_grade_csum[1] = -cond_logodds_loh_ovc/2;
  cond_logodds_loh_grade_csum[2] = cond_logodds_loh_ovc/2;
  for (i in 3:n_grade) {
    cond_logodds_loh_grade_csum[i] = cond_logodds_loh_grade_csum[i-1] + cond_logodds_loh_grade[i-2] * sd_cond_logodds_loh_grade + mean_cond_logodds_loh_grade;
  }
  
  //construct positive scale vector
  vector[n_grade] cond_prob_loh = inv_logit(cond_logodds_loh_base + cond_logodds_loh_grade_csum);
  
  //find each of the mixture probabilities
  array[n_grade, 5] real mix_prob;
  for (i in 1:n_grade) {
    mix_prob[i,1] = prob_AB[i]; //probability in allelic balance
    mix_prob[i,2] = ((1-prob_AB[i]) * cond_prob_loh[i]) / 2; //probability not in AB but in LoH (tail 1)
    mix_prob[i,3] = mix_prob[i,2]; //probability not in AB but in LoH (tail 2)
    mix_prob[i,4] = ((1-prob_AB[i]) * (1-cond_prob_loh[i])) / 2; //probability not in AB but in slab (tail 1)
    mix_prob[i,5] = mix_prob[i,4]; //probability not in AB but in slab (tail 2)
  }
  array[n_grade, 5] real log_mix_prob = log(mix_prob);
  
}

model {
  
  //////////////////////////// 
  // priors and hyperpriors //
  //////////////////////////// 
  
  mu_base ~ normal((bounds[1] + bounds[2])/2, slab_width/6);
  mu_ovc ~ normal(0, slab_width/4);
  mean_mu_grade ~ normal(0, slab_width/4);
  mu_grade ~ std_normal();
  sd_mu_grade ~ normal(0,0.2);
  
  sigma_base ~ std_normal();
  sigma_ovc ~ std_normal();
  mean_sigma_grade ~ std_normal();
  sigma_grade ~ std_normal();
  sd_sigma_grade ~ std_normal();
  
  logodds_AB_base ~ std_normal();
  logodds_AB_ovc ~ normal(0, 0.5);
  mean_logodds_AB_grade ~ normal(0, 0.5);
  logodds_AB_grade ~ std_normal();
  sd_logodds_AB_grade ~ normal(0, 0.5);
  
  cond_logodds_loh_base ~ std_normal();
  cond_logodds_loh_ovc ~ normal(0, 0.5);
  mean_cond_logodds_loh_grade ~ normal(0, 0.5);
  cond_logodds_loh_grade ~ std_normal();
  sd_cond_logodds_loh_grade ~ normal(0, 0.5);
  
  sigma_theta_loh ~ normal(0, 0.5);
  sigma_theta_AB ~ normal(0, 0.5);
  
  theta_AB ~ double_exponential(0.5, sigma_theta_AB);
  theta_slab ~ normal(mu[grade_k], sigma[grade_k]);
  theta_loh ~ exponential(1/sigma_theta_loh);
  
  ////////////////////////
  // Mixture Likelihood //
  ////////////////////////
  
  //average binomial probabilities across mixture
  for (i in 1:n) {
    array[5] real lp;
    int i_grade = grade_k[indiv_j[i]];
    lp[1] = log_mix_prob[i_grade,1] + binomial_lpmf(count[i] | total[i], theta_AB[indiv_j[i]]);
    lp[2] = log_mix_prob[i_grade,2] + binomial_lpmf(count[i] | total[i], theta_loh[indiv_j[i]]);
    lp[3] = log_mix_prob[i_grade,3] + binomial_lpmf(count[i] | total[i], 1-theta_loh[indiv_j[i]]);
    lp[4] = log_mix_prob[i_grade,4] + binomial_lpmf(count[i] | total[i], theta_slab[indiv_j[i]]);
    lp[5] = log_mix_prob[i_grade,5] + binomial_lpmf(count[i] | total[i], 1-theta_slab[indiv_j[i]]);
    target += log_sum_exp(lp);
  }
  
}
