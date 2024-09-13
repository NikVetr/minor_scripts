
data{
    int p;
    matrix[p, 2] x_err;
    matrix<lower=0>[p, 2] sd_x_err;
}
parameters{
    array[p] vector[2] x;
    vector<lower=0>[2] sigma_x;
    cholesky_factor_corr[2] L;
}
model{
    sigma_x ~ normal(0,2);
    x_err[,1] ~ normal(x[,1], sd_x_err[,1]);
    x_err[,2] ~ normal(x[,2], sd_x_err[,2]);
    L ~ lkj_corr_cholesky(1);
    x ~ multi_normal_cholesky(rep_vector(0, 2), diag_pre_multiply(sigma_x, L));
}
generated quantities{
    real r = multiply_lower_tri_self_transpose(L)[1,2];
}
