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
  real<lower = -sqrt(1-ss_obs), upper = sqrt(1-ss_obs)> x_unobs;
  //vector<lower=-1, upper=1>[m] z;    // Unobserved sample from unit hypersphere
  
}

transformed parameters {
  vector[m] x;
  x[1:(m-1)] = x_obs;
  x[m] = x_unobs;
  //transform away from sqrt(y) * z
  //vector[m] z_sq = (x / sqrt(y))^2;
  vector[m] z_sq = x^2 / y; // hmm, this does not observe the simplex constraint
  //what if I used stick breaking to specify a beta on the last element of z_sq?
}

model {
  // Prior for y from Beta distribution with m/2 and alpha
  y ~ beta(m / 2.0, alpha);

  // Implicit prior for z to be from unit hypersphere
  z_sq ~ dirichlet(rep_vector(0.5, m)); 
  
  // need jacobian adjustment for z_sq
  // implying x_obs ~ sqrt(y) * z[1:(m-1)]; // this is not how this is done
  target += m * log(2) + sum(log(abs(x))) - m * log(y);

}

generated quantities {

  // very last element (pth, pth) of the cholesky factor
  real xp = sqrt(1 - y);
  
  // correlation from cholesky factor product
  real r = sum(rpm1 .* x);
}
