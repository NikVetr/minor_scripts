functions {
  //-- log-pdf for p independent 2-D observations with row-specific diag errors
  real bv_norm_me_lpdf(matrix x,              // p × 2 observed values
                        vector mu,             // length-2 mean
                        vector sigma,          // length-2 SD of latent
                        real   r,              // correlation of latent
                        vector sd1,            // length-p meas-error SD col-1
                        vector sd2) {          // length-p meas-error SD col-2
    int    p = rows(x);
    vector[p] a  = square(sigma[1]) + square(sd1);     // Σ₁₁ + σ²ₑ,₁
    vector[p] b  = square(sigma[2]) + square(sd2);     // Σ₂₂ + σ²ₑ,₂
    real        c  = r * sigma[1] * sigma[2];          // Σ₁₂ = Σ₂₁

    vector[p] det = a .* b - square(c);                // |Vᵢ|
    vector[p] dx1 = x[,1] - mu[1];
    vector[p] dx2 = x[,2] - mu[2];
    vector[p] quad = ( b .* square(dx1)
                     + a .* square(dx2)
                     - 2 * c .* (dx1 .* dx2) ) ./ det; // (x-μ)'V⁻¹(x-μ)

    return -p * log(2 * pi())              // constant part (d = 2)
           -0.5 * sum(log(det) + quad);    // variable part
  }
}

data {
  int<lower=1>  p;
  matrix[p,2]   x_err;
  matrix<lower=0>[p,2] sd_x_err;
}
parameters {
  vector[2]        mu;
  vector<lower=0>[2] sigma;
  real<lower=-1,upper=1> r;
}
model {
  // priors
  mu    ~ normal(0,5);
  sigma ~ normal(0,10);
  r     ~ uniform(-1,1);

  // likelihood (single line, no loops, no extra parameters)
  target += bv_norm_me_lpdf(x_err | mu, sigma, r,
                             sd_x_err[,1], sd_x_err[,2]);
}
