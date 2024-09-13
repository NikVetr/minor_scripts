
data {
  int<lower=1> N;
  int<lower=1> D;
  array[N] vector[D] y;
}

transformed data {
  int n_combn = choose(D,2);
  array[n_combn, 2] int combm;
  int i = 1;
  for(j in 1:D){
    for(k in (j+1):D){
        combm[i,1] = j;
        combm[i,2] = k;
        i += 1;
    }
  }
}

parameters {
  vector<lower=-1, upper=1>[n_combn] r;
}

model {
  
  //implicit prior on r
  (r + 1) / 2 ~ beta(D / 2.0, D / 2.0);
  
  //likelihood calculation
  for(ic in 1:n_combn){
    
    //get data
    array[N] vector[2] y_pair;
    y_pair[,1] = y[,combm[ic,1]];
    y_pair[,2] = y[,combm[ic,2]];
    
    //construct correlation matrix
    matrix[2,2] R_pair;
    R_pair[1,2] = r[ic];
    R_pair[2,1] = r[ic];
    R_pair[1,1] =  1;
    R_pair[2,2] =  1;
    
    //evaluate bivariate likelihood
    y_pair ~ multi_normal(rep_vector(0, 2), R_pair);
  }
  
}

generated quantities {
  matrix[D, D] R; //can't be a corr_matrix because PSD not enforced
  for(ic in 1:n_combn){
    R[combm[ic,1],combm[ic,2]] = r[ic];
    R[combm[ic,2],combm[ic,1]] = r[ic];
  }
}

