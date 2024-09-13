
functions {
  int sum2d(array[,] int a) {
    int s = 0;
    for (i in 1:size(a)) {
      s += sum(a[i]);
    }
    return s;
  }
}

data {
  int<lower=1> N;
  int<lower=1> D;
  array[N, D] int<lower=0, upper=1> y;
}

transformed data {
  int<lower=0> N_pos;
  array[sum2d(y)] int<lower=1, upper=N> n_pos;
  array[size(n_pos)] int<lower=1, upper=D> d_pos;
  int<lower=0> N_neg;
  array[(N * D) - size(n_pos)] int<lower=1, upper=N> n_neg;
  array[size(n_neg)] int<lower=1, upper=D> d_neg;

  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i;
    int j;
    i = 1;
    j = 1;
    for (n in 1:N) {
      for (d in 1:D) {
        if (y[n, d] == 1) {
          n_pos[i] = n;
          d_pos[i] = d;
          i += 1;
        } else {
          n_neg[j] = n;
          d_neg[j] = d;
          j += 1;
        }
      }
    }
  }
}

parameters {
  vector[D] mu;
  cholesky_factor_corr[D] L;
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
}

transformed parameters {
  array[N] vector[D] z;
  for (n in 1:N_pos) {
    z[n_pos[n], d_pos[n]] = z_pos[n];
  }
  for (n in 1:N_neg) {
    z[n_neg[n], d_neg[n]] = z_neg[n];
  }
}

model {
  L ~ lkj_corr_cholesky(1);
  to_vector(mu) ~ normal(0, 5);
  z ~ multi_normal_cholesky(mu, L);
}

generated quantities {
  corr_matrix[D] R;
  R = multiply_lower_tri_self_transpose(L);
}

