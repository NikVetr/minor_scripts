data {
  int<lower=1> n0; //number of initial individuals
  int<lower=1> R; //number of filtration rounds
  array[R] int<lower=1> n; //number of individuals left after each filtration round
  int<lower=1> G; //number of group categories
  array[G] int<lower=1> g; //number of groups in each group category
  array[n0, G] int<lower=1> gi; //individual group index in each group category
  array[choose(g, 2), 2] int<lower=1> gpk; //order of group pairings
  array[n0] int fri; //number of rounds progressed by each individual
}

transformed data{
  int<lower=1> npairs = choose(g, 2);
  int<lower=1> maxg = max(g); //for ragged matrices (to have buffer, can fix to 0)
}

parameters {
  
  //define main effects
  real B_intercept; 
  vector[R] B_round; 
  vector[n0] B_indiv; 
  array[G, max(g)] B_group;
  
  //define interaction effects
  array[npairs, maxg, maxg] B_group_x_group;
  array[npairs, maxg, R] B_group_x_round;
  array[npairs, maxg, maxg, R] B_group_x_group_x_round;
  
  //define scale
  real<lower=0> B_round_sd;
  real<lower=0> B_indiv_sd;
  real<lower=0> B_group_sd;
  real<lower=0> B_group_x_group_sd;
  real<lower=0> B_group_x_round_sd;
  real<lower=0> B_group_x_group_x_round_sd;

}

model {
  //priors
  
  //likelihood
  
  //iterate through each individuals, and all the rounds they individually progressed through
  //constructing their probability of success until the final round where they met failure
  for(i in 1:n0){
    vector[fri[i]] r_liab = rep_vector(B_indiv[i], fri[i]) + B_intercept; //initialize with their individual effect
    r_liab += B_round[1:fri[i]]; //add the main effects for the round
    
    //add the group effects for all rounds
    for(k in 1:G){
      r_liab += B_group[k, gi[i,k]];
    }
    
    //add the group interaction effects for all rounds
    for(k in 1:npairs){
      r_liab += B_group_x_group[k, gi[i,gpk[k,1]], gi[i,gpk[k,2]]];
    }
    
    //add the round-wise effects for group and group interactions
    for(j in 1:fri[i]){
      for(k in 1:G){
        
      }
      
      //add the group interaction effects for all rounds
      for(k in 1:npairs){
          
      }
    }
  }
}

