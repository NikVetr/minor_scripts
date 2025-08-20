data {
  int<lower=2> k;                  // size in 1st dim   (rows)
  int<lower=2> p;                  // size in 2nd dim   (columns)
  int<lower=2> l;                  // size in 3rd dim   (slices)

  matrix[k, k-1] Q_k;              // zero‑sum basis, dim‑1
  matrix[p, p-1] Q_p;              // zero‑sum basis, dim‑2
  matrix[l, l-1] Q_l;              // zero‑sum basis, dim‑3

  int<lower=1,upper=4> prior;      // 1=N, 2=Laplace, 3=t(ν=2), 4=Cauchy
}

transformed data {
  // variance‑matching factor: 1 / sqrt(prod_d (1-1/n_d))
  real sd_scale_factor =
        inv_sqrt((1 - 1.0 / k) * (1 - 1.0 / p) * (1 - 1.0 / l));
  int<lower=1> kp1 = (k - 1) * (p - 1) * (l - 1);   // # free params
}

parameters {
  vector[kp1] theta_raw_flat;      // vectorised   (k‑1)(p‑1)(l‑1)

  real<lower=0> tau;               // global scale (needed by 2–4)
}

transformed parameters {

  //------------------------------------------------
  // 4.1  reshape the flat vector → 3‑D array
  //------------------------------------------------
  array[k-1, p-1, l-1] theta_raw;
  {
    int pos = 1;
    for (i in 1:(k-1))
      for (j in 1:(p-1))
        for (m in 1:(l-1)) {
          theta_raw[i, j, m] = theta_raw_flat[pos];
          pos += 1;
        }
  }

  //------------------------------------------------
  // 4.2  apply the latent scale  (prior‑specific)
  //------------------------------------------------
  array[k-1, p-1, l-1] theta_scaled;

  if      (prior == 1)  theta_scaled = theta_raw;                 // Normal
  else if (prior == 2)  theta_scaled = sqrt(tau) * theta_raw;     // Laplace
  else if (prior == 3)  theta_scaled = theta_raw / sqrt(tau);     // Student‑t
  else if (prior == 4)  theta_scaled = theta_raw / sqrt(tau);     // Cauchy

  //------------------------------------------------
  // 4.3  project into zero‑sum sub‑space
  //      Θ = s * (Q_k  θ  Q_pᵀ)   along dim‑1 & dim‑2,
  //      then multiply along dim‑3 by Q_lᵀ
  //------------------------------------------------
  array[k, p-1, l-1] tmp1;           // after dim‑1 projection
  for (j in 1:(p-1))
    for (m in 1:(l-1))
      tmp1[:, j, m] = Q_k * theta_scaled[:, j, m];

  array[k, p, l-1]  tmp2;            // after dim‑2 projection
  {
    vector[p] colbuf;
    for (i in 1:k)
      for (m in 1:(l-1)) {
        for (j in 1:(p-1))
          colbuf[j] = tmp1[i, j, m];
        tmp2[i, :, m] = Q_p * head(colbuf, p - 1);
      }
  }

  array[k, p, l] theta;              // after dim‑3 projection
  {
    vector[l] slice_buf;
    for (i in 1:k)
      for (j in 1:p) {
        for (m in 1:(l-1))
          slice_buf[m] = tmp2[i, j, m];
        theta[i, j, :] = sd_scale_factor * (Q_l * head(slice_buf, l - 1));
      }
  }
}

model {

  //---------------- base normals ------------------
  to_vector(theta_raw_flat) ~ std_normal();

  //---------------- latent scales -----------------
  if (prior == 2) {                 // Laplace(0,1)
      tau ~ exponential(0.5);       // rate = 1/(2 b²), b = 1
  } else if (prior == 3) {          // Student‑t(ν = 2)
      tau ~ gamma(1, 1);
  } else if (prior == 4) {          // Cauchy (t₁)
      tau ~ gamma(0.5, 0.5);
  }
}

