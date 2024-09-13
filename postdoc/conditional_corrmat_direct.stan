data {
  int<lower=1> m;
  real<lower=0> eta;
  matrix<lower=-1,upper=1>[m,m] R_obs;         // eta shape parameter
  vector<lower=-1,upper=1>[2] r_bounds;         // Observed components of x
}

parameters {
  real<lower=r_bounds[1], upper=r_bounds[2]> r;
}

transformed parameters {
  matrix<lower=-1,upper=1>[m,m] R;
  for(i in 1:m){
    
    for(j in 1:m){
      
      int R_observed = 0;
      
      if(i == j){
        R[i,j] = 1;
        R_observed += 1;
      }
      
      if(i == (m-1) && j == m){
        R[i,j] = r;
        R_observed += 1;
      }
      
      if(j == (m-1) && i == m){
        R[i,j] = r;
        R_observed += 1;
      }
      
      if(R_observed == 0) {
        R[i,j] = R_obs[i,j];
      }
      
    }
  }
}

model {
  R ~ lkj_corr(eta);
}
