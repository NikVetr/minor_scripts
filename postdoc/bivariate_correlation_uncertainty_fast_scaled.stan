data {
  int<lower=1>  p;
  matrix[p,2]   x_err;
  matrix<lower=0>[p,2] sd_x_err;
}
parameters {
  array[p] vector[2] z;        // standard normals
  vector<lower=0>[2] sigma;
  real<lower=-1,upper=1> r;
}
transformed parameters {
  matrix[2,2] Sigma;                 // latent covariance
  Sigma[1,1] = square(sigma[1]);
  Sigma[2,2] = square(sigma[2]);
  Sigma[1,2] = r * sigma[1] * sigma[2];
  Sigma[2,1] = Sigma[1,2];
  array[p] vector[2] x;        // the true values
  {
    matrix[2,2] InvS = inverse(Sigma);
    for (i in 1:p) {
      matrix[2,2] Dinv = diag_matrix(to_vector(inv(square(sd_x_err[i]))));
      matrix[2,2] V    = inverse(InvS + Dinv);     // conditional var
      vector[2]  m     = V * (Dinv * x_err[i]');       // conditional mean
      x[i] = m + cholesky_decompose(V) * z[i];
    }
  }
}
model {
  z[,1]  ~ std_normal();
  z[,2]  ~ std_normal();
  sigma         ~ normal(0,2);
  r             ~ uniform(-1,1);
  // *no* likelihood term — it is absorbed in the transformation Jacobian
}
