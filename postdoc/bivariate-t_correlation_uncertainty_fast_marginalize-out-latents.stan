functions {
  /**
   * Log-pdf for p independent 2-D observations with
   *   - latent N(mu, Sigma)             (integrated out),
   *   - Student-t measurement error     (handled via gamma scales).
   *
   * x        : p × 2 matrix of observed values
   * mu       : length-2 mean vector
   * sigma    : latent SDs (length-2, >0)
   * r        : latent correlation (|r|<1)
   * sd1, sd2 : length-p measurement-error scales
   * lam1,lam2: length-p gamma precision draws (λ)
   */
  real bv_t_me_lpdf(matrix x,
                    vector mu,
                    vector sigma,
                    real   r,
                    vector sd1,
                    vector sd2,
                    vector lam1,
                    vector lam2) {
    int p = rows(x);

    // row-specific diagonal elements of V_i = Σ + diag(s^2 / λ)
    vector[p] a = square(sigma[1]) + square(sd1) ./ lam1;
    vector[p] b = square(sigma[2]) + square(sd2) ./ lam2;
    real       c = r * sigma[1] * sigma[2];

    vector[p] det  = a .* b - square(c);            // |V_i|
    vector[p] dx1  = x[,1] - mu[1];
    vector[p] dx2  = x[,2] - mu[2];
    vector[p] quad = ( b .* square(dx1)
                     + a .* square(dx2)
                     - 2 * c .* (dx1 .* dx2) ) ./ det;

    return -p * log(2 * pi())
           -0.5 * sum(log(det) + quad);
  }
}

data {
  int<lower=1>           p;
  matrix[p,2]            x_err;
  matrix<lower=0>[p,2]   sd_x_err;
  array[p,2] int<lower=1> df;          // d.f. estimated with x_err
}

parameters {
  // 4 global latent-covariance parameters
  vector[2]              mu;
  vector<lower=0>[2]     sigma;
  real<lower=-1,upper=1> r;

  // 2 × p gamma precision draws for the Student-t error
  vector<lower=0>[p]     lam1;
  vector<lower=0>[p]     lam2;
}

model {
  // -- gamma priors implementing the Student-t mixture
  lam1 ~ gamma(to_vector(df[,1]) ./ 2, to_vector(df[,1]) ./ 2);
  lam2 ~ gamma(to_vector(df[,2]) ./ 2, to_vector(df[,2]) ./ 2);

  // usual weakly-informative priors for the latent field
  mu    ~ normal(0, 5);
  sigma ~ normal(0,10);
  r     ~ uniform(-1,1);

  // likelihood (no explicit loop)
  target += bv_t_me_lpdf(x_err | mu, sigma, r,
                         sd_x_err[,1], sd_x_err[,2], lam1, lam2);
}
