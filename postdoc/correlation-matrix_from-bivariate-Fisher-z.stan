data{
  int p;
  int np;
  array[np,2] int ip;
  vector[np] z;
  vector[np] sigma_z;
}

parameters{
  //corr_matrix[p] R;
  cholesky_factor_corr[p] L;
}

transformed parameters {
  corr_matrix[p] R = tcrossprod(L);
  vector[np] Rz;
  for(i in 1:np){
    real r = R[ip[i,1], ip[i,2]];
    Rz[i] = 0.5 * log( (1 + r) / (1 - r) );
  }
}

model{
  z ~ normal(Rz, sigma_z);
  //R ~ lkj_corr(1);
  L ~ lkj_corr_cholesky(1);
}

