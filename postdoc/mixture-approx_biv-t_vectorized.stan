functions {
  // element-wise log-sum-exp of two vectors
  vector lse2(vector a, vector b) {
    int N = rows(a);
    vector[N] out;
    for (n in 1:N) {
      real m = fmax(a[n], b[n]);                    // scalar max
      out[n] = m + log(exp(a[n] - m) + exp(b[n] - m));
    }
    return out;
  }
}

data {
  int<lower=1>  p;
  matrix[p,2]   x_err;
  matrix<lower=0>[p,2] sd_x_err;
  int<lower=1>  K;                  // K = 3

  // mixture tables  (p rows × K components)
  matrix[p,K]   w1;     // weights for column 1   (sum to 1 row–wise)
  matrix[p,K]   w2;     // weights for column 2
  matrix[p,K]   s2_1;   // extra variances for col-1   ( = 1/λ_k )
  matrix[p,K]   s2_2;   // extra variances for col-2
}

parameters {
  vector[2]              mu;
  vector<lower=0>[2]     sigma;
  real<lower=-1,upper=1> r;
}

model {
  //---------------- priors
  mu    ~ normal(0,5);
  sigma ~ normal(0,10);
  r     ~ uniform(-1,1);

  //---------------- data-dependent helpers
  vector[p] dx1 = x_err[,1] - mu[1];
  vector[p] dx2 = x_err[,2] - mu[2];
  real       c  = r * sigma[1] * sigma[2];

  //---------------- first mixture component (k = 1)
  vector[p] a = square(sd_x_err[,1]) + s2_1[,1] + square(sigma[1]);
  vector[p] b = square(sd_x_err[,2]) + s2_2[,1] + square(sigma[2]);
  vector[p] det  = a .* b - square(c);
  vector[p] quad = ( b .* square(dx1)
                   + a .* square(dx2)
                   - 2 * c .* (dx1 .* dx2) ) ./ det;

  // joint log-weight  log(w1_i1 * w2_i1)
  vector[p] log_mix = log( w1[,1] .* w2[,1] ) - 0.5 * (log(det) + quad);

  //---------------- remaining K-1 components
  for (k in 2:K) {
    a    = square(sd_x_err[,1]) + s2_1[,k] + square(sigma[1]);
    b    = square(sd_x_err[,2]) + s2_2[,k] + square(sigma[2]);
    det  = a .* b - square(c);
    quad = ( b .* square(dx1)
           + a .* square(dx2)
           - 2 * c .* (dx1 .* dx2) ) ./ det;

    vector[p] log_comp = log( w1[,k] .* w2[,k] )
                         - 0.5 * (log(det) + quad);
    log_mix = lse2(log_mix, log_comp);     // element-wise log-sum-exp
  }

  target += -p * log(2*pi()) + sum(log_mix);
}
