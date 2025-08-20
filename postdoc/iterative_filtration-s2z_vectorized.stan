functions {
  /**
   * project a (K-1) vector to a sum-to-zero K-vector
   * (unit variance once multiplied by its own sd parameter)
   */
  vector proj_vec(matrix Q, vector raw) {
    int K = rows(Q);
    return (Q * raw) * inv_sqrt(1 - 1.0 / K);
  }

  /**
   * fast 1-axis projection of a flat tensor
   *  z      – full vector to overwrite in place
   *  z_core – indices of the (K-1)×N core
   *  z_full – indices of the K×N full matrix
   *  Qk     – K × (K-1) zero-sum basis for the first axis
   *  K1m    – K-1
   *  N      – N  (∏ other dims)
   */
  void project_axis1(inout vector z,
                     array[] int z_core,
                     array[] int z_full,
                     matrix Qk,
                     int K1m,
                     int N) {
    matrix[K1m, N] core = to_matrix(z[z_core], K1m, N);
    matrix[K1m+1, N] full = Qk * core;
    z[z_full] = to_vector(full);
  }
}

data {
  // ---------- global sizes ----------
  int<lower=1> R;                     // rounds
  array[3] int<lower=2> g;            // groups per category
  int<lower=1> N;                     // total Bernoulli trials (rows in long table)

  // ---------- look-up columns ----------
  array[N]   int<lower=0,upper=1>   y;
  array[N]   int<lower=1>           r_id;        // still used for B_round
  array[N]   int<lower=1> idx_g1r;
  array[N]   int<lower=1> idx_g2r;
  array[N]   int<lower=1> idx_g3r;
  array[N]   int<lower=1> idx_12;
  array[N]   int<lower=1> idx_13;
  array[N]   int<lower=1> idx_23;
  array[N]   int<lower=1> idx_12r;
  array[N]   int<lower=1> idx_13r;
  array[N]   int<lower=1> idx_23r;

  // ---------- zero-sum tensor indices (pre-computed) ----------
  //  pairwise g × g
  array[] int z12;   array[] int z12c;   int ncol12;
  array[] int z13;   array[] int z13c;   int ncol13;
  array[] int z23;   array[] int z23c;   int ncol23;
  //  pairwise × R   (stored as g1×g2×R, etc.)
  array[] int z12r;  array[] int z12rc;  int ncol12r;   // ncol = g2*R
  array[] int z13r;  array[] int z13rc;  int ncol13r;
  array[] int z23r;  array[] int z23rc;  int ncol23r;

  // ---------- round effect indices ----------
  int<lower=1> P_round;           // = R
  array[P_round] int<lower=1> z_round;    // 1 … R (identity, but keeps code symmetric)

  // ---------- Q bases ----------
  matrix[R,   R-1]     Q_R;
  matrix[g[1], g[1]-1] Q1;
  matrix[g[2], g[2]-1] Q2;
  matrix[g[3], g[3]-1] Q3;
}

parameters {
  real           B_intercept;

  // ----- round main effect -------
  vector[R-1]    raw_round;      real<lower=0> sd_round;

  // ----- main groups ------------
  vector[g[1]-1] raw_g1;         real<lower=0> sd_g1;
  vector[g[2]-1] raw_g2;         real<lower=0> sd_g2;
  vector[g[3]-1] raw_g3;         real<lower=0> sd_g3;

  // ----- group × round ----------
  vector[(g[1]-1)*(R-1)] raw_g1r;
  vector[(g[2]-1)*(R-1)] raw_g2r;
  vector[(g[3]-1)*(R-1)] raw_g3r;
  real<lower=0> sd_gx_round;

  // ----- group × group ----------
  vector[(g[1]-1)*(g[2]-1)] raw_12;
  vector[(g[1]-1)*(g[3]-1)] raw_13;
  vector[(g[2]-1)*(g[3]-1)] raw_23;
  real<lower=0> sd_gxg;

  // ----- group × group × round ---
  vector[(g[1]-1)*(g[2]-1)*(R-1)] raw_12r;
  vector[(g[1]-1)*(g[3]-1)*(R-1)] raw_13r;
  vector[(g[2]-1)*(g[3]-1)*(R-1)] raw_23r;
  real<lower=0> sd_gxg_x_round;
}

transformed parameters {
  // ---------- main round ----------
  vector[R] B_round = sd_round * proj_vec(Q_R, raw_round);

  // ---------- main groups ----------
  vector[g[1]] B_g1 = sd_g1 * proj_vec(Q1, raw_g1);
  vector[g[2]] B_g2 = sd_g2 * proj_vec(Q2, raw_g2);
  vector[g[3]] B_g3 = sd_g3 * proj_vec(Q3, raw_g3);

  // ---------- group × round (flattened gᵢ × R matrices) ----------
  vector[g[1]*R] B_g1r = rep_vector(0, g[1]*R);
  {
    int K1m = g[1]-1;
    project_axis1(B_g1r, z1c, z1, Q1, K1m, ncol1);   // reuse z1 indices
    B_g1r *= sd_gx_round;
  }

  vector[g[2]*R] B_g2r = rep_vector(0, g[2]*R);
  {
    int K1m = g[2]-1;
    project_axis1(B_g2r, z2c, z2, Q2, K1m, ncol2);
    B_g2r *= sd_gx_round;
  }

  vector[g[3]*R] B_g3r = rep_vector(0, g[3]*R);
  {
    int K1m = g[3]-1;
    project_axis1(B_g3r, z3c, z3, Q3, K1m, ncol3);
    B_g3r *= sd_gx_round;
  }

  // ---------- group × group ----------
  vector[g[1]*g[2]] B12  = rep_vector(0, g[1]*g[2]);
  {
    int K1m = g[1]-1;
    project_axis1(B12, z12c, z12, Q1, K1m, ncol12);
    B12  *= sd_gxg;
  }

  vector[g[1]*g[3]] B13  = rep_vector(0, g[1]*g[3]);
  project_axis1(B13, z13c, z13, Q1, g[1]-1, ncol13);
  B13 *= sd_gxg;

  vector[g[2]*g[3]] B23  = rep_vector(0, g[2]*g[3]);
  project_axis1(B23, z23c, z23, Q2, g[2]-1, ncol23);
  B23 *= sd_gxg;

  // ---------- group × group × round ----------
  vector[g[1]*g[2]*R] B12r = rep_vector(0, g[1]*g[2]*R);
  {
    int K1m = g[1]-1;
    project_axis1(B12r, z12rc, z12r, Q1, K1m, ncol12r);
    B12r *= sd_gxg_x_round;
  }

  vector[g[1]*g[3]*R] B13r = rep_vector(0, g[1]*g[3]*R);
  project_axis1(B13r, z13rc, z13r, Q1, g[1]-1, ncol13r);
  B13r *= sd_gxg_x_round;

  vector[g[2]*g[3]*R] B23r = rep_vector(0, g[2]*g[3]*R);
  project_axis1(B23r, z23rc, z23r, Q2, g[2]-1, ncol23r);
  B23r *= sd_gxg_x_round;
}

model {
  // ---- raw priors (vectorised) ----
  raw_round  ~ std_normal();
  raw_g1     ~ std_normal();
  raw_g2     ~ std_normal();
  raw_g3     ~ std_normal();
  raw_g1r    ~ std_normal();
  raw_g2r    ~ std_normal();
  raw_g3r    ~ std_normal();
  raw_12     ~ std_normal();
  raw_13     ~ std_normal();
  raw_23     ~ std_normal();
  raw_12r    ~ std_normal();
  raw_13r    ~ std_normal();
  raw_23r    ~ std_normal();

  // ---- scale priors --------------
  sd_round         ~ exponential(1);
  sd_g1 | sd_g2 | sd_g3  ~ exponential(1);
  sd_gx_round      ~ exponential(1);
  sd_gxg           ~ exponential(1);
  sd_gxg_x_round   ~ exponential(1);

  B_intercept      ~ std_normal();

  // ---------- vectorised likelihood ----------
  vector[N] eta = rep_vector(B_intercept, N);

  // main round
  eta += B_round[r_id];

  // main groups
  eta += B_g1[g1_id] + B_g2[g2_id] + B_g3[g3_id];

  // group × round
  eta += B_g1r[idx_g1r] + B_g2r[idx_g2r] + B_g3r[idx_g3r];

  // group × group
  eta += B12[idx_12] + B13[idx_13] + B23[idx_23];

  // group × group × round
  eta += B12r[idx_12r] + B13r[idx_13r] + B23r[idx_23r];

  y ~ bernoulli_logit(eta);
}

generated quantities {
  real lp__ = target();            // for debugging / loo
}
