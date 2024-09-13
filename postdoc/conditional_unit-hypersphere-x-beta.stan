data {
  int<lower=1> m;            
  real<lower=0> eta;         // eta shape parameter
  vector<lower=-1,upper=1>[m-1] x_obs;         // Observed components of x
  vector<lower=-1,upper=1>[m] rpm1;            // previous column of the cholesky factor (to compute r in genquan block)
}

transformed data {
  real alpha = eta + (m - 1) / 2.0 - 0.5 * m;  // Calculate alpha
  real ss_obs = dot_self(x_obs);  // sub of squares of the observed components of x
}

parameters {
  real<lower=1-ss_obs, upper=1> y;       // The scaling factor from Beta(m/2, alpha)
  //real<lower = -sqrt(1-ss_obs), upper = sqrt(1-ss_obs)> x_unobs;
  real<lower = 0, upper = sqrt(1-ss_obs)> x_unobs; //only focus on positive region
  
  //the last hypersphere element is sqrt(beta(1/2, (k-1)/2)) distributed, with reflection across the origin
  
}

transformed parameters {
  //transform away from sqrt(y) * z
  real zm = x_unobs / sqrt(y); 
  real zm_sqrt = sqrt(zm);
  
}

model {
  // Prior for y from Beta distribution with m/2 and alpha
  y ~ beta(m / 2.0, alpha);

  // Implicit prior for z to be from unit hypersphere
  zm_sqrt ~ beta(0.5, (m - 1.0) / 2.0);
  
  // need jacobian adjustment for zm_sqrt
  // implying x_obs ~ sqrt(y) * z[1:(m-1)]; // (this is not how this is done)
  target += - 0.5 * log(x_unobs) - 0.25 * log(y);


}

generated quantities {
  
  vector[m] x;
  x[1:(m-1)] = x_obs;
  x[m] = x_unobs * ((binomial_rng(1, 0.5) * 2.0) - 1.0);
  
  // very last element (pth, pth) of the cholesky factor
  real xp = sqrt(1 - y);
  
  // correlation from cholesky factor product
  real r = sum(rpm1 .* x);
}
