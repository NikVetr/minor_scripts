functions {

  real guess_theta(real theta_q, real w1, real w2, real w3, real mu, real sigma, real lambda) {
    real x;
    
    if (theta_q < w2) {
      // Left truncated exponential component
      real u = theta_q / w2;  // Rescaled quantile within this component
      real Z_E = 1 - exp(-lambda * 1);  // Normalizing constant for truncation
      x = - (1 / lambda) * log1m(u * Z_E);
    } else if (theta_q < (w2 + w1)) {
      // Middle truncated normal component
      real u = (theta_q - w2) / w1;  // Rescaled quantile within this component
      real Phi_L = normal_cdf(0 | mu, sigma);
      real Phi_U = normal_cdf(1 | mu, sigma);
      real inv_Phi_arg = u * (Phi_U - Phi_L) + Phi_L;
      // Ensure inv_Phi_arg is within (0, 1)
      inv_Phi_arg = fmin(1 - 1e-10, fmax(1e-10, inv_Phi_arg));
      x = mu + sigma * inv_Phi(inv_Phi_arg);
      
    } else {
      // Right reflected truncated exponential component
      real u_rel = (theta_q - w2 - w1) / w3;  // Rescaled quantile within this component
      real exp_neg_lambda = exp(-lambda);
      real Z_E = 1 - exp_neg_lambda;  // Normalizing constant for truncation
      real exp_arg = u_rel * Z_E + exp_neg_lambda;
      // Ensure exp_arg is positive
      exp_arg = fmax(1e-10, exp_arg);
      x = 1 + (1 / lambda) * log(exp_arg);
    }
    
      return x;
  }
  
  vector mixture_qf(vector x, vector theta_vec, array[] real x_r, array[] int x_i) {
    // Extract parameters from theta_vec
    real mean_normal = theta_vec[1];
    real stddev_normal = theta_vec[2];
    real rate_exponential = theta_vec[3];
    real weight_normal = theta_vec[4];
    real weight_exponential = theta_vec[5];
    real weight_ref_truncexp = theta_vec[6];
    real target_quantile = theta_vec[7];

    real current_theta = inv_logit(x[1]);

    // Compute truncated normal CDF
    real truncnorm_total_mass = normal_cdf(1 | mean_normal, stddev_normal) - 
                           normal_cdf(0 | mean_normal, stddev_normal);
    real truncnorm_mass = (normal_cdf(current_theta | mean_normal, stddev_normal) - 
                           normal_cdf(0 | mean_normal, stddev_normal)) / truncnorm_total_mass;

    // Compute truncated exponential CDF
    real truncexp_total_mass = 1 - exp(-rate_exponential);
    real truncexp_mass = (1 - exp(-rate_exponential * current_theta)) / 
                                     truncexp_total_mass;

    // Compute reflected truncated exponential CDF
    real ref_truncexp_mass = 1 - ((1 - exp(-rate_exponential * (1 - current_theta))) / truncexp_total_mass);

    // Combine into mixture CDF
    real mixture_mass = weight_normal * truncnorm_mass +
                     weight_exponential * truncexp_mass +
                     weight_ref_truncexp * ref_truncexp_mass;

    // Initialize result vector
    vector[1] result;
    result[1] = mixture_mass - target_quantile;
    
    return result;
}

}

data {
  real mu;     // Mean of the normal distribution
  real sigma;  // Standard deviation of the normal distribution
  real lambda; // Rate of the exponential distributions
  real w1;     // Weight of the normal component
  real w2;     // Weight of the exponential component
  real w3;     // Weight of the reflected exponential component
  int<lower=1, upper=2> solver_index; //1 uses newton, 2 uses powell
}

transformed data {
  array[0] real x_r;
  array[0] int x_i;
  real rel_tol = 1E-12;
  real f_tol = 1E-3;
  real scaling_step = 1E-3;
  int max_steps = 5000;
}

parameters {
  real<lower=0, upper=1> theta_q;  // Uniform(0, 1) quantile
}

transformed parameters {
  real theta;
  vector[7] theta_vec;
  vector[1] theta_guess;
  theta_guess[1] = logit(fmin(fmax(guess_theta(theta_q, w1, w2, w3, mu, sigma, lambda), 1E-6), 1-1E-6));

  theta_vec[1] = mu;
  theta_vec[2] = sigma;
  theta_vec[3] = lambda;
  theta_vec[4] = w1;
  theta_vec[5] = w2;
  theta_vec[6] = w3;
  theta_vec[7] = theta_q;
  
  if(solver_index == 1){
    theta = inv_logit(solve_newton_tol(mixture_qf, theta_guess, 
                                       scaling_step, f_tol, max_steps,
                                       theta_vec, x_r, x_i)[1]);
  } else {
    theta = inv_logit(solve_powell_tol(mixture_qf, theta_guess, 
                                       rel_tol, f_tol, max_steps,
                                       theta_vec, x_r, x_i)[1]);
  }
  
  // print("theta_q = ", theta_q);
  // print("theta_guess = ", inv_logit(theta_guess));
  // print("theta = ", theta);
  // print("-------");
  
}

model {
  // Implicitly using Uniform(0, 1) prior for u
}