functions {
  int count_elem(int[] test, int elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
      if(test[i] == elem)
        count = count + 1;
    return(count);
  }
  
  int[] which_elem(int[] test, int elem) {
    int res[count_elem(test, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(test))
      if(test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
}

data {
    int<lower=1> n;
    int<lower=1> n_prot;
    int<lower=1> n_group;
    int<lower=1> n_time;
    int<lower=1> n_pat;
    int<lower=0> censor_threshold;
    //int<lower=0> count[n];
    real log_count[n];
    int<lower=0, upper=1> uncens[n];
    int<lower=1, upper=n_group> group[n];
    int<lower=1, upper=n_time> time[n];
    int<lower=1, upper=n_prot> prot[n];
    int<lower=1, upper=n_pat> pat[n];
}

transformed data {
    real log_censor_threshold = log(censor_threshold);
    
    int<lower=0> n_cens = count_elem(uncens, 0);
    int cens_idx[n_cens] = which_elem(uncens, 0);

    int<lower=0> n_uncens = count_elem(uncens, 1);
    int uncens_idx[n_uncens] = which_elem(uncens, 1);
    real uncens_log_count[n_uncens] = log_count[uncens_idx];
}

parameters {
    //censored observations
    real<upper=log_censor_threshold> cens_log_count[n_cens];
    
    //model parameters
    real B_raw[n_time-1, n_group, n_prot];
    real<lower=0> B_sd[n_time-1, n_group];
    real B_mean[n_time-1, n_group];
    real<lower=0> sigma[n_time-1];
}

transformed parameters{
    //fill in both obs and censored (sampled) obs
    real all_log_counts[n_time, n_group, n_prot, n_pat];
    real diff_log_counts[n_time-1, n_group, n_prot, n_pat];
    real B[n_time-1, n_group, n_prot];


    //fill in all_log_counts
    // for(t in 1:n_time){
    //   all_log_counts[t,,,] = rep_array(0.0, n_group, n_prot, n_pat); //initialize array
    // }
    
    for(i in 1:n_uncens){
      all_log_counts[time[uncens_idx[i]], 
                     group[uncens_idx[i]], 
                     prot[uncens_idx[i]], 
                     pat[uncens_idx[i]]] = uncens_log_count[i];
    }

    for(i in 1:n_cens){
      all_log_counts[time[cens_idx[i]], 
                     group[cens_idx[i]], 
                     prot[cens_idx[i]], 
                     pat[cens_idx[i]]] = cens_log_count[i];
    }

    //eval difference between the log counts
    // for(t in 1:(n_time-1)){
    //   diff_log_counts[t,,,] = rep_array(0.0, n_group, n_prot, n_pat); //initialize array
    // }
    
    for(t in 1:(n_time-1)){
      for(g in 1:n_group){
          diff_log_counts[t,g,,] = to_array_2d(to_matrix(all_log_counts[t+1,g,,]) - to_matrix(all_log_counts[t,g,,]));
      }
    }

    //center the parameters
    for(t in 1:(n_time-1)){
      for(g in 1:n_group){
        for(p in 1:n_prot){
          B[t,g,p] = B_mean[t,g] + B_raw[t,g,p] * B_sd[t,g];
        }
      }
    }

}
model {
    //(hyper)priors for coef
    for(t in 1:(n_time-1)){
      B_mean[t,] ~ normal(0,10);
      B_sd[t,] ~ normal(0,2);
      for(g in 1:n_group){
        B_raw[t,g,] ~ std_normal();  
      }
    }
    sigma ~ normal(0,2);
    
    //likelihood
    for(t in 1:(n_time-1)){
      for(g in 1:n_group){
        for(p in 1:n_prot){
          diff_log_counts[t,g,p,] ~ normal(B[t,g,p], to_vector(rep_array(sigma[t],n_pat)));
        }
      }
    }
}