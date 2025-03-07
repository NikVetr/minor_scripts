data {
  //for first blob (eg left eye)
  vector[2] blob_1_mu;
  vector<lower=0>[2] blob_1_sd;
  
  //for second blob (eg right eye)
  vector[2] blob_2_mu;
  vector<lower=0>[2] blob_2_sd;
  
  //for banana (eg mouth)
  vector[2] banana_mu;
  vector<lower=0>[2] banana_sd;
  
  //probabilities for all components
  vector[3] probs;
}
transformed data {
  real log_p_banana = log(probs[1]);
  real log_p_blob_1   = log(probs[2]);
  real log_p_blob_2   = log(probs[3]);

}
parameters {
  real x;
  real y;
}
model {
  //component densities
  real lp_banana = normal_lpdf(x | banana_mu[1], banana_sd[1])
                 + normal_lpdf(y | banana_mu[2] + x^2, banana_sd[2]);
  real lp_blob_1  = normal_lpdf(x | blob_1_mu[1], blob_1_sd[1])
                 + normal_lpdf(y | blob_1_mu[2], blob_1_sd[2]);
  real lp_blob_2  = normal_lpdf(x | blob_2_mu[1], blob_2_sd[1])
                 + normal_lpdf(y | blob_2_mu[2], blob_2_sd[2]);

  //mixture density:
  vector[3] lp;
  lp[1] = log_p_banana + lp_banana;
  lp[2] = log_p_blob_1 + lp_blob_1;
  lp[3] = log_p_blob_2 + lp_blob_2;
  target += log_sum_exp(lp);
}
