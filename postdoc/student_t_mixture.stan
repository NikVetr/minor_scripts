functions {
  // log‑density of a K‑component normal mixture with shared mean
  real mix_normal_lpdf(real y, vector pi, real mu, vector sigma) {
    int K = num_elements(pi);
    vector[K] log_comp;
    for (k in 1:K)
      log_comp[k] = log(pi[k]) + normal_lpdf(y | mu, sigma[k]);
    return log_sum_exp(log_comp);
  }
}

data {
  int<lower=1> K;            // mixture components
  real<lower=0> df;          // Student‑t d.o.f.
  real loc;                  // shared location (mean)
  real<lower=0> scale;       // scale of target t
  vector[K] lambda;          // quadrature nodes (Gamma precisions)
  vector[K] pi;              // quadrature weights (sum to 1)
}

transformed data {
  vector[K] sigma_comp;
  for (k in 1:K)
    sigma_comp[k] = scale / sqrt(lambda[k]);
}

parameters {
  real x;                    // draw from exact t prior
  real y;                    // draw from mixture prior
}

model {
  // exact Student‑t prior for x
  target += student_t_lpdf(x | df, loc, scale);
  
  // finite Gaussian‑mixture prior for y (latent k is marginalised)
  target += mix_normal_lpdf(y | pi, loc, sigma_comp);
}
