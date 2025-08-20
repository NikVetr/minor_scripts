data {
  int<lower=2> k;                       // length of theta
  matrix[k, k-1] Q_k;                   // zero‑sum orthonormal basis
  int<lower=1> prior; 
}

transformed data {
  real sd_scale_factor = inv_sqrt(1 - 1.0 / k);
}

parameters {
  vector[k-1] theta_raw;                  // latent vector (always size k‑1)

  // scale for the mixture priors
  real<lower=0> tau;
}

transformed parameters {
  
  vector[k-1] theta_scaled;
  
  // rescale the iid standard normal vector to some other distribution
  if      (prior == 1)  theta_scaled = theta_raw;                   // Normal
  else if (prior == 2)  theta_scaled = sqrt(tau) * theta_raw;       // Laplace
  else if (prior == 3)  theta_scaled = theta_raw / sqrt(tau);       // Student‑t
  else if (prior == 4)  theta_scaled = theta_raw / sqrt(tau);       // Cauchy (0, ν=1)
  else if (prior == 5)  theta_scaled = theta_raw * sqrt(tau);       // NEG(0, 2, 1)
  
  // project to orthonormal manifold and rescale
  vector[k] theta = sd_scale_factor * (Q_k * theta_scaled);
  
}

model {

  // unconstrained variable always gets an iid standard normal prior
  theta_raw ~ std_normal();
  
  // priors for rescaling auxiliary variables

  // give tau a std_normal prior to avoid divergent iterations (never gets used)
  if (prior == 1) {
    tau ~ std_normal();
    
  // double_exponential(0,1)
  } else if (prior == 2) {
      real rate = 1.0;
      tau ~ exponential(1 / (2 * rate^2));
  
  // student_t(2, 0, 1)
  } else if (prior == 3) {
      real nu = 2.0;
      tau ~ gamma(nu / 2, nu / 2);
  
  // cauchy(0,1)
  } else if (prior == 4) {
      real nu = 1.0;
      tau ~ gamma(nu / 2, nu / 2);
  
    // NEG(0, 2, 1)
  } else if (prior == 5) {
      real lambda = 1;
      real alpha = 2;
      tau ~ pareto_type_2(0, lambda, alpha);
  }
}
