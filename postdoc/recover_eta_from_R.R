library(cmdstanr)
library(posterior)
library(caret)
library(MASS)

#sample full dimension correlation matrix
p <- 10
eta <- 3
R <- rlkj(p)
uR <- chol(R)

#build Stan model
stan_loc <- "~/scripts/minor_scripts/postdoc/conditional_unit-hypersphere-x-beta.stan"
stan_loc <- "~/scripts/minor_scripts/postdoc/conditional_corrmat_direct.stan"
stan_program <- paste0(readLines(stan_loc), collapse = "\n")
mod <- cmdstan_model(stan_loc)

#build data object
dat <- list(m = p - 1,
            eta = 1,
            x_obs = uR[1:(p-2),p],
            rpm1 = uR[1:(p-1),p-1])

dat <- list(m = p,
            eta = 1,
            R_obs = R,
            r_bounds = unlist(coef_bounds(R)))


fit <- mod$sample(chains = 4, iter_sampling = 5E3, iter_warmup = 5E3, data = dat, 
                  parallel_chains = 4, adapt_delta = 0.9, max_treedepth = 10, 
                  refresh = 100, init = 0.1)
# summ <- fit$summary("missing_r")
samps <- as.data.frame(as_draws_df(fit$draws("missing_r")))
hist(samps$missing_r, breaks = 20)