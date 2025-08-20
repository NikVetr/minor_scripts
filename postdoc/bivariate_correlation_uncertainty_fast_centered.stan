
data{
    int p;
    matrix[p, 2] x_err;
    matrix<lower=0>[p, 2] sd_x_err;
}
parameters{
    vector[p] x1_raw;
    vector[p] x2_raw;
    vector<lower=0>[2] sigma_x;
    real<lower=-1, upper=1> r;
}

transformed parameters{
  matrix[p, 2] x;
  x[,1] = x1_raw;
  x[,2] = r .* x1_raw + sqrt(1-square(r)) .* x2_raw;
}

model{
    //prior
    sigma_x ~ normal(0,2);
    x1_raw ~ normal(0, sigma_x[1]);
    x2_raw ~ normal(0, sigma_x[2]);
    
    //likelihood
    x_err[,1] ~ normal(x[,1], sd_x_err[,1]);
    x_err[,2] ~ normal(x[,2], sd_x_err[,2]);
}
