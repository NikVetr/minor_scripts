data {
  int<lower=1>  p;
  matrix[p,2]   x_err;
  matrix<lower=0>[p,2] sd_x_err;
}
parameters {
  vector[2] mu;          // latent mean
  vector<lower=0>[2] sigma;          // latent s.d.
  real<lower=-1,upper=1> r;          // latent correlation
}
transformed parameters {
  matrix[2,2] Sigma;                 // latent covariance
  Sigma[1,1] = square(sigma[1]);
  Sigma[2,2] = square(sigma[2]);
  Sigma[1,2] = r * sigma[1] * sigma[2];
  Sigma[2,1] = Sigma[1,2];
}
model {
  mu ~ normal(0,5);
  sigma ~ normal(0, 10);
  r     ~ uniform(-1, 1);
  for (i in 1:p) {
    // Add measurement error for this row
    matrix[2,2] V = Sigma;
    V[1,1] += square(sd_x_err[i,1]);
    V[2,2] += square(sd_x_err[i,2]);
    x_err[i]   ~ multi_normal(mu, V);
  }
}
