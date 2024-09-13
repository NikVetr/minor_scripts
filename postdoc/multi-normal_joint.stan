
data {
  int<lower=1> N;
  int<lower=1> D;
  array[N] vector[D] y;
}

parameters {
  cholesky_factor_corr[D] L;
}

model {
  L ~ lkj_corr_cholesky(1);
  y ~ multi_normal_cholesky(rep_vector(0,D), L);
}

generated quantities {
  corr_matrix[D] R;
  R = multiply_lower_tri_self_transpose(L);
}

