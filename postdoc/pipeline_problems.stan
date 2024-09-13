data {
  int<lower=1> n0; //number of initial individuals
  int<lower=1> R; //number of filtration rounds
  array[R+1] int<lower=1> n; //number of individuals left after each filtration round
  int<lower=1> G; //number of group categories
  array[G] int<lower=1> g; //number of groups in each group category
  array[n0, G] int<lower=1> gi; //individual group index in each group category
  array[choose(G, 2), 2] int<lower=1> gpk; //order of group pairings
  array[n0] int fri; //number of rounds progressed by each individual
}

transformed data{
  int<lower=1> npairs = choose(G, 2);
  int<lower=1> maxg = max(g); //for ragged matrices (to have buffer, can fix to 0)
}

parameters {
  
  //define main effects
  real B_intercept; 
  vector[R] B_round; 
  vector[n0] B_indiv; 
  array[G, max(g)] real B_group;
  
  //define interaction effects
  array[npairs, maxg, maxg] real B_group_x_group;
  array[G, maxg, R] real B_group_x_round;
  array[npairs, maxg, maxg, R] real B_group_x_group_x_round;
  
  //define scale
  real<lower=0> B_round_sd;
  real<lower=0> B_indiv_sd;
  real<lower=0> B_group_sd;
  real<lower=0> B_group_x_group_sd;
  real<lower=0> B_group_x_round_sd;
  real<lower=0> B_group_x_group_x_round_sd;
  
  //TODO : define all these to have further multilevel effects, eg group-wise scale

}

model {
  /////////////////////////////////
  //////////// priors ////////////
  ////////////////////////////////
  
  // main effects
  B_intercept ~ std_normal();
  B_round ~ std_normal(); 
  B_indiv ~ std_normal(); 
  for(k in 1:G){
    B_group[k,] ~ std_normal();  
  }
  
  // interactions
  for(k in 1:npairs){
    for(j in 1:maxg){
      B_group_x_group[k,j,] ~ std_normal();
      for(l in 1:maxg){
        B_group_x_group_x_round[k,j,l,] ~ std_normal();    
      }
    }
  }
  
  for(k in 1:G){
    for(j in 1:maxg){
      B_group_x_round[k,j,] ~ std_normal();
    }
  }
  
  // scale
  B_round_sd ~ std_normal();
  B_indiv_sd ~ std_normal();
  B_group_sd ~ std_normal();
  B_group_x_group_sd ~ std_normal();
  B_group_x_round_sd ~ std_normal();
  B_group_x_group_x_round_sd ~ std_normal();
  
  ////////////////////////////////////
  //////////// likelihood ////////////
  ////////////////////////////////////
  
  //iterate through each individual, and for all the rounds they individually progressed through
  //compose the probability of their success in each round until the final round where they met failure
  for(i in 1:n0){
    
    //how far did the person make it?
    array[2] int nrs;
    nrs[1] = fri[i];
    nrs[2] = R;
    int nr = min(nrs);
    array[nr] int track_record = rep_array(1, nr); 
    if(fri[i] != (R+1)){
      track_record[nr] = 0;
    }
    
    vector[nr] r_liab = rep_vector(B_indiv[i], nr) * B_indiv_sd + B_intercept; //initialize with their individual effect
    r_liab += B_round[1:nr] * B_round_sd * 100; //add the main effects for the round
    
    //add the group effects for all rounds
    for(k in 1:G){
      r_liab += rep_vector(B_group[k, gi[i,k]] * B_group_sd, nr);
    }
    
    //add the group interaction effects for all rounds
    for(k in 1:npairs){
      r_liab += rep_vector(B_group_x_group[k, gi[i,gpk[k,1]], gi[i,gpk[k,2]]] * B_group_x_group_sd, nr);
    }
    
    //add the round-wise effects for group and group x group interactions
    for(j in 1:nr){
      for(k in 1:G){
        r_liab[j] += B_group_x_round[k, gi[i,k], j] * B_group_x_round_sd;
      }
      
      //add the group interaction effects for each round
      for(k in 1:npairs){
          r_liab[j] += B_group_x_group_x_round[k, gi[i,gpk[k,1]], gi[i,gpk[k,2]], j] * B_group_x_group_x_round_sd;
      }
    }
    //r_liab now contains the liabilities for that individual moving on through their rounds
    
    //compute likelihood of individual's history of successes and failure
    track_record ~ binomial(rep_array(1, nr), inv_logit(r_liab));
    
  }
}

