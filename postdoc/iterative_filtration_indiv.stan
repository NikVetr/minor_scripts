data {
  int<lower=1> n0; //number of initial individuals
  int<lower=1> R; //number of filtration rounds
  array[R+1] int<lower=1> n; //number of individuals left after each filtration round
  int<lower=1> G; //number of group categories
  array[G] int<lower=1> g; //number of groups in each group category
  array[n0, G] int<lower=1> gi; //individual group index in each group category
  array[choose(G, 2), 2] int<lower=1> gpk; //order of group pairings
  array[n0] int<lower=1,upper=R+1> fri; //number of rounds progressed by each individual
}

transformed data{
  array[n0] int obs_r; //how many rounds of selection did we observe for this person?
  for(i in 1:n0){
    obs_r[i] = min(fri[i], R);
  }
}

parameters {
  // intercept, round main, individual random effect
  real B_intercept;
  vector[R]   z_round;
  vector[n0]  z_indiv;
  real<lower=0> sd_round;
  real<lower=0> sd_indiv;

  // main group effects
  vector[g[1]] z_g1;
  vector[g[2]] z_g2;
  vector[g[3]] z_g3;
  real<lower=0> sd_g1;
  real<lower=0> sd_g2;
  real<lower=0> sd_g3;

  // group × round (one matrix per category)
  matrix[g[1], R] z_g1_x_round;
  matrix[g[2], R] z_g2_x_round;
  matrix[g[3], R] z_g3_x_round;
  real<lower=0> sd_gx_round;

  // group × group
  matrix[g[1], g[2]] z_g1x2;
  matrix[g[1], g[3]] z_g1x3;
  matrix[g[2], g[3]] z_g2x3;
  real<lower=0> sd_gxg;

  // group × group × round
  array[R] matrix[g[1], g[2]] z_g1x2_x_round;
  array[R] matrix[g[1], g[3]] z_g1x3_x_round;
  array[R] matrix[g[2], g[3]] z_g2x3_x_round;
  real<lower=0> sd_gxg_x_round;
}

transformed parameters {
  vector[R]  B_round = sd_round * z_round;
  vector[n0] B_indiv = sd_indiv * z_indiv;

  vector[g[1]] B_g1 = sd_g1 * z_g1;
  vector[g[2]] B_g2 = sd_g2 * z_g2;
  vector[g[3]] B_g3 = sd_g3 * z_g3;

  matrix[g[1], R] B_g1_x_round = sd_gx_round * z_g1_x_round;
  matrix[g[2], R] B_g2_x_round = sd_gx_round * z_g2_x_round;
  matrix[g[3], R] B_g3_x_round = sd_gx_round * z_g3_x_round;

  matrix[g[1], g[2]] B_g1x2 = sd_gxg * z_g1x2;
  matrix[g[1], g[3]] B_g1x3 = sd_gxg * z_g1x3;
  matrix[g[2], g[3]] B_g2x3 = sd_gxg * z_g2x3;

  array[R] matrix[g[1], g[2]] B_g1x2_x_round;
  array[R] matrix[g[1], g[3]] B_g1x3_x_round;
  array[R] matrix[g[2], g[3]] B_g2x3_x_round;
  for (r in 1:R) {
 B_g1x2_x_round[r] = sd_gxg_x_round * z_g1x2_x_round[r];
 B_g1x3_x_round[r] = sd_gxg_x_round * z_g1x3_x_round[r];
 B_g2x3_x_round[r] = sd_gxg_x_round * z_g2x3_x_round[r];
  }
}

model {
  // priors
  B_intercept ~ std_normal();
  z_round  ~ std_normal();
  z_indiv  ~ std_normal();
  z_g1 ~ std_normal();
  z_g2 ~ std_normal();
  z_g3 ~ std_normal();
  to_vector(z_g1_x_round) ~ std_normal();
  to_vector(z_g2_x_round) ~ std_normal();
  to_vector(z_g3_x_round) ~ std_normal();
  to_vector(z_g1x2) ~ std_normal();
  to_vector(z_g1x3) ~ std_normal();
  to_vector(z_g2x3) ~ std_normal();
  for (r in 1:R) {
   to_vector(z_g1x2_x_round[r]) ~ std_normal();
   to_vector(z_g1x3_x_round[r]) ~ std_normal();
   to_vector(z_g2x3_x_round[r]) ~ std_normal();
  }

  sd_round  ~ std_normal();
  sd_indiv  ~ std_normal();
  sd_g1 ~ std_normal();
  sd_g2 ~ std_normal();
  sd_g3 ~ std_normal();
  sd_gx_round  ~ std_normal();
  sd_gxg ~ std_normal();
  sd_gxg_x_round ~ std_normal();

  // likelihood
 for (i in 1:n0) {
   
   // get information on the observed rounds for this person
   int nr = obs_r[i]; // rounds observed for this person
   
   // also y, their observed success history: 1’s until failure, 0 at end if they
   // did not make it all the way
   array[nr] int y = rep_array(1, nr);
   if(fri[i] <= R){
    y[nr] = 0; // they failed the last round
   }
   
   // initialize container for progress liability
   vector[nr] eta = B_intercept + B_indiv[i] + B_round[1:nr];
  
   // add main group effects
   eta += rep_vector(B_g1[gi[i,1]], nr);
   eta += rep_vector(B_g2[gi[i,2]], nr);
   eta += rep_vector(B_g3[gi[i,3]], nr);
  
   // add group × round
   for (r in 1:nr) {
    eta[r] += B_g1_x_round[gi[i,1], r];
    eta[r] += B_g2_x_round[gi[i,2], r];
    eta[r] += B_g3_x_round[gi[i,3], r];
   }
  
   // add group × group
   eta += rep_vector(B_g1x2[gi[i,1], gi[i,2]], nr);
   eta += rep_vector(B_g1x3[gi[i,1], gi[i,3]], nr);
   eta += rep_vector(B_g2x3[gi[i,2], gi[i,3]], nr);
  
   // add group × group × round
   for (r in 1:nr) {
    eta[r] += B_g1x2_x_round[r][gi[i,1], gi[i,2]];
    eta[r] += B_g1x3_x_round[r][gi[i,1], gi[i,3]];
    eta[r] += B_g2x3_x_round[r][gi[i,2], gi[i,3]];
   }
    
   y ~ bernoulli_logit(eta);
  }
}
