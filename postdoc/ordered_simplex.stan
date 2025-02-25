functions {

  // Beta(1,b) CDF:  F(x) = 1 - (1-x)^b
  // x in [0,1], b>0
  real beta1b_cumdf(real x, real b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    return 1 - pow(1 - x, b);
  }

  // Inverse CDF for Beta(1,b):  F^{-1}(q) = 1 - (1 - q)^(1/b)
  real beta1b_icdf(real q, real b) {
    if (q <= 0) return 0;
    if (q >= 1) return 1;
    return 1 - pow(1 - q, 1/b);
  }

}

data {
  int<lower=2> K;  // dimension
}

parameters {
  // K-1 uniform(0,1) to drive the stepwise truncated Beta draws
  vector<lower=0, upper=1>[K-1] u;
  simplex[K] y_raw;
}

transformed parameters {
  vector[K] y = sort_desc(y_raw);
  vector[K] x;
  real remaining = 1.0;  // leftover "mass"
  real x_prev = 1.0;     // treat x_0=1 for bounding x_1
  
  vector[K-1] lost_mass;
  for (k in 1:(K-1)) {
    // 1) Lower bound: lb = remaining/(K - k + 1)
    real lb = remaining / (K - k + 1);

    // 2) Upper bound: ub = min(x_{k-1}, remaining), with x_{0} = 1
    real ub = fmin(x_prev, remaining);

    // 3) Convert to fractional bounds in [0,1]
    //    fraction = x_k / remaining
    real lbFrac = lb / remaining;
    real ubFrac = ub / remaining;

    // Beta(1,b) with b = (K-k)
    real b = K - k;

    // 4) Evaluate Beta(1,b) cdf at lbFrac, ubFrac
    real cdf_lb = beta1b_cumdf(lbFrac, b);
    real cdf_ub = beta1b_cumdf(ubFrac, b);
    lost_mass[k] = cdf_ub - cdf_lb;
    
    // 5) sample quantile uniformly in [cdf_lb, cdf_ub]
    real quantile = cdf_lb + u[k] * (cdf_ub - cdf_lb);

    // 6) inverse cdf => fraction in [lbFrac, ubFrac]
    real frac_k = beta1b_icdf(quantile, b);

    // 7) scale by leftover => x_k
    x[k] = frac_k * remaining;

    // 8) update leftover and "previous" for next iteration
    remaining -= x[k];
    x_prev = x[k];
  }

  // final coordinate
  x[K] = remaining;
}

model {
  // No likelihood increments => purely a forward transform from u -> x.
  target += -log(1-lost_mass);
}
