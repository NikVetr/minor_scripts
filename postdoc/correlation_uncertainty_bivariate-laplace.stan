
data{
    int n;
    matrix[n, 2] x_err;
    matrix<lower=0>[n, 2] sd_x_err;
}
parameters{
    array[n] vector[2] x;
    vector<lower=0>[n] w; // auxiliary variable to get bivariate exponential
    vector<lower=0>[2] sigma_x;
    cholesky_factor_corr[2] L;
}
model{
    sigma_x ~ normal(0,2);
    w ~ exponential(1);
    x_err[,1] ~ normal(to_vector(x[,1]) .* sqrt(w) , sd_x_err[,1]);
    x_err[,2] ~ normal(to_vector(x[,2]) .* sqrt(w), sd_x_err[,2]);
    L ~ lkj_corr_cholesky(1);
    x ~ multi_normal_cholesky(rep_vector(0, 2), diag_pre_multiply(sigma_x, L));
}
generated quantities{
    real r = multiply_lower_tri_self_transpose(L)[1,2];
}
