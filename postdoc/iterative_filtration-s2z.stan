functions {
 // project a (K-1)-vector onto sum-to-zero K-vector, rescale to unit variance
 vector proj_vec(matrix Q, vector raw) {
   int K = dims(Q)[1];
   real scale = inv_sqrt(1.0 - 1.0 / K); // ≡ sqrt((K-1)/K)
   return scale * (Q * raw);
 }

 // project a (K-1)×(P-1) matrix onto K×P double-centred, unit-var matrix
 matrix proj_mat(matrix Qk, matrix raw, matrix Qp) {
   int K = dims(Qk)[1];
   int P = dims(Qp)[1];
   real scale = inv_sqrt( (1.0 - 1.0/K) * (1.0 - 1.0/P) );
   return scale * (Qk * raw * Qp');
 }
}

data {
 int<lower=1> n0; //number of initial individuals
 int<lower=1> R; //number of filtration rounds
 array[R+1] int<lower=1> n; //number of individuals left after each filtration round
 int<lower=1> G; //number of group categories
 array[G] int<lower=1> g; //number of groups in each group category
 array[n0, G] int<lower=1> gi; //individual group index in each group category
 array[choose(G, 2), 2] int<lower=1> gpk; //order of group pairings
 array[n0] int<lower=1,upper=R+1> fri; //number of rounds progressed by each individual
 matrix[R, R-1] Q_R;
 matrix[g[1], g[1]-1] Q_g1;
 matrix[g[2], g[2]-1] Q_g2;
 matrix[g[3], g[3]-1] Q_g3;
}

transformed data{
 array[n0] int obs_r; //how many rounds of selection did we observe for this person?
 for(i in 1:n0){
 obs_r[i] = min(fri[i], R);
 }
}

parameters {
 real B_intercept;

 vector[R-1] z_round_raw;
 real<lower=0> sd_round;

 // main groups -
 vector[g[1]-1] z_g1_raw;
 vector[g[2]-1] z_g2_raw;
 vector[g[3]-1] z_g3_raw;
 real<lower=0> sd_g1;
 real<lower=0> sd_g2;
 real<lower=0> sd_g3;

 // group × round -
 matrix[g[1]-1, R-1] z_g1_x_round_raw;
 matrix[g[2]-1, R-1] z_g2_x_round_raw;
 matrix[g[3]-1, R-1] z_g3_x_round_raw;
 real<lower=0> sd_gx_round;

 // group × group -
 matrix[g[1]-1, g[2]-1] z_g1x2_raw;
 matrix[g[1]-1, g[3]-1] z_g1x3_raw;
 matrix[g[2]-1, g[3]-1] z_g2x3_raw;
 real<lower=0> sd_gxg;

 // group × group × round (one raw matrix per round) 
 array[R-1] matrix[g[1]-1, g[2]-1] z_g1x2_x_round_raw;
 array[R-1] matrix[g[1]-1, g[3]-1] z_g1x3_x_round_raw;
 array[R-1] matrix[g[2]-1, g[3]-1] z_g2x3_x_round_raw;
 real<lower=0> sd_gxg_x_round;
}

transformed parameters {

 // project vectors -
 vector[R] B_round = sd_round * proj_vec(Q_R, z_round_raw);
 vector[g[1]] B_g1 = sd_g1 * proj_vec(Q_g1, z_g1_raw);
 vector[g[2]] B_g2 = sd_g2 * proj_vec(Q_g2, z_g2_raw);
 vector[g[3]] B_g3 = sd_g3 * proj_vec(Q_g3, z_g3_raw);

 // project matrices 
 matrix[g[1],R] B_g1_x_round = sd_gx_round * proj_mat(Q_g1, z_g1_x_round_raw, Q_R);
 matrix[g[2],R] B_g2_x_round = sd_gx_round * proj_mat(Q_g2, z_g2_x_round_raw, Q_R);
 matrix[g[3],R] B_g3_x_round = sd_gx_round * proj_mat(Q_g3, z_g3_x_round_raw, Q_R);

 matrix[g[1],g[2]] B_g1x2 = sd_gxg * proj_mat(Q_g1, z_g1x2_raw, Q_g2);
 matrix[g[1],g[3]] B_g1x3 = sd_gxg * proj_mat(Q_g1, z_g1x3_raw, Q_g3);
 matrix[g[2],g[3]] B_g2x3 = sd_gxg * proj_mat(Q_g2, z_g2x3_raw, Q_g3);

 // project group × group × round -
 array[R] matrix[g[1],g[2]] B_g1x2_x_round;
 array[R] matrix[g[1],g[3]] B_g1x3_x_round;
 array[R] matrix[g[2],g[3]] B_g2x3_x_round;
 
 // center each pairwise group slice in tensor
 array[R-1] matrix[g[1],g[2]] tmp12;
 array[R-1] matrix[g[1],g[3]] tmp13;
 array[R-1] matrix[g[2],g[3]] tmp23;

 for (r in 1:(R-1)) {
  tmp12[r] = proj_mat(Q_g1, z_g1x2_x_round_raw[r], Q_g2);
  tmp13[r] = proj_mat(Q_g1, z_g1x3_x_round_raw[r], Q_g3);
  tmp23[r] = proj_mat(Q_g2, z_g2x3_x_round_raw[r], Q_g3);
 }

 //now center each fiber across rounds
 for (i in 1:g[1]) {
  for (j in 1:g[2]) {
    vector[R-1] raw_vec = rep_vector(0, R-1);
    for (r in 1:(R-1))  raw_vec[r] = tmp12[r][i,j];
    vector[R] centred = proj_vec(Q_R, raw_vec) * sd_gxg_x_round;
    for (r in 1:R)  B_g1x2_x_round[r][i,j] = centred[r];
  }
 }
 for (i in 1:g[1]) {
  for (j in 1:g[3]) {
    vector[R-1] raw_vec = rep_vector(0, R-1);
    for (r in 1:(R-1))  raw_vec[r] = tmp13[r][i,j];
    vector[R] centred = proj_vec(Q_R, raw_vec) * sd_gxg_x_round;
    for (r in 1:R)  B_g1x3_x_round[r][i,j] = centred[r];
  }
 }
 for (i in 1:g[2]) {
  for (j in 1:g[3]) {
    vector[R-1] raw_vec = rep_vector(0, R-1);
    for (r in 1:(R-1))  raw_vec[r] = tmp23[r][i,j];
    vector[R] centred = proj_vec(Q_R, raw_vec) * sd_gxg_x_round;
    for (r in 1:R)  B_g2x3_x_round[r][i,j] = centred[r];
  }
 }

}

model {
 // raw priors
 z_round_raw ~ std_normal();
 z_g1_raw ~ std_normal();
 z_g2_raw ~ std_normal();
 z_g3_raw ~ std_normal();
 to_vector(z_g1_x_round_raw) ~ std_normal();
 to_vector(z_g2_x_round_raw) ~ std_normal();
 to_vector(z_g3_x_round_raw) ~ std_normal();
 to_vector(z_g1x2_raw) ~ std_normal();
 to_vector(z_g1x3_raw) ~ std_normal();
 to_vector(z_g2x3_raw) ~ std_normal();
 for (r in 1:(R-1)) {
   to_vector(z_g1x2_x_round_raw[r]) ~ std_normal();
   to_vector(z_g1x3_x_round_raw[r]) ~ std_normal();
   to_vector(z_g2x3_x_round_raw[r]) ~ std_normal();
 }

 // scales
 sd_round ~ exponential(1);
 sd_g1 ~ exponential(1);
 sd_g2 ~ exponential(1);
 sd_g3 ~ exponential(1);
 sd_gx_round ~ exponential(1);
 sd_gxg ~ exponential(1);
 sd_gxg_x_round ~ exponential(1);

 // grand intercept
 B_intercept ~ std_normal();

 // likelihood 
 for (i in 1:n0) {
   int nr = obs_r[i];
   vector[nr] eta = B_intercept + B_round[1:nr];

   eta += rep_vector(B_g1[gi[i,1]], nr);
   eta += rep_vector(B_g2[gi[i,2]], nr);
   eta += rep_vector(B_g3[gi[i,3]], nr);

   for (r in 1:nr) {
     eta[r] += B_g1_x_round[gi[i,1], r];
     eta[r] += B_g2_x_round[gi[i,2], r];
     eta[r] += B_g3_x_round[gi[i,3], r];
   }

   eta += rep_vector(B_g1x2[gi[i,1], gi[i,2]], nr);
   eta += rep_vector(B_g1x3[gi[i,1], gi[i,3]], nr);
   eta += rep_vector(B_g2x3[gi[i,2], gi[i,3]], nr);

   for (r in 1:nr) {
     eta[r] += B_g1x2_x_round[r][gi[i,1], gi[i,2]];
     eta[r] += B_g1x3_x_round[r][gi[i,1], gi[i,3]];
     eta[r] += B_g2x3_x_round[r][gi[i,2], gi[i,3]];
   }

   array[nr] int y = rep_array(1, nr);
   if (fri[i] <= R) y[nr] = 0;
   y ~ bernoulli_logit(eta);
 }
}
