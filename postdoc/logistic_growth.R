library(rethinking); set_ulam_cmdstan(F)
library(cmdstanr)
library(rstan)
library(splines)

#### simulate data ####

#true model specification, P = K / (1+((K-L0)/L0)e^(-rt))
set.seed(123)
P <- function(K, r, t, L0){
  return(K / (1 + ((K-L0)/L0)*exp(-r*t) ))
}
K <- 2500 #asymptote of growth process
t0 = 0; t1 = 64 #start and end time of the growth process
t <- t0:(t1*4)/4
r = log(K) / (t1/2) #set r such that flex point is halfway through time range
L0 <- 1 #starting state
rates <- P(K, r, t, L0) #expected rates at each timestep
n_obs <- 128 #number of observations to draw
prop_to_sample <- 1
sampleTimes <- sort(runif(round(n_obs * prop_to_sample), range(t)[1], range(t)[2] * prop_to_sample)) #uniformly distributed sampling times
sampleExpRates <- P(K, r, sampleTimes, L0) #sampled expected rates
theta <- 8 #scale parameter of gamma distribution
sampleRates <- rgamma(round(n_obs * prop_to_sample), shape = sampleExpRates / theta, scale = theta) #sampled rates at each timepoint
sampleObs <- rpois(round(n_obs * prop_to_sample), sampleRates) #poisson-distributed counts at each timepoint
d <- list(time = sampleTimes, count = sampleObs) #put data into list
plot(t, rates, type = "l", ylim = c(0,5000))
points(sampleTimes, sampleObs)

#### describe animation ####

# phase 0 -- have arnold and frank talk "your method of inverse probability won't avail you here" using speech bubbles and xkcd font popping in one at a time, 
# efron labrada can give fisher the bootstrap, franco newton an apple
# have them comment throughout, e.g. "incredible" pun on credible intervals, "so uninformed" on flat priors, "A MAP by any other name" for penalized likelihood, etc.
# bayes can pop in all "Hi folks! Today we'll be exploring how soon into a non-linear temporal process it's end behavior can be well understood. We'll be using a linear combination of logistic functions as our non-linear model of choice, featuring both growth and decay components (draw these) whose sum determine (sum them) the system's overall behavior. 
# These yield an expected rate, which we use to parameterize a negative-binomial distribution to realize rates by which poisson-distributed counts may be drawn (do the lightning bolt thing). We'll fit the true, perfectly specified model to these simulated data one observation at a time to show you the true power of Applied Bayesian Sta-"
# then Fisher rises from the bottom of the frame saying "Not so fast, BARNIE! Your paltry method of inverse probability will not avail you here! Truly, it is only the Method of Maximum Likelihood that provides a principled means of estimating model parameters."
# then he'll pause and say "Wait, how do I express the logistic MLE in closed form again? Little help here boys!" and newton and efron and turing will appear saying "You have my legacy!" "And my method!" "And MY bootstraps!"
#we simulate sampling times in that most natural of ways, by drawing from a rescaled dirichlet(5,5,5...)


# phase 1a -- fullscreen plot, visualize curve growing, visualize exponentially distributed sampling events 
# visualize gamma distribution moving, visualize poisson distribution forming, shoot arrows from poisson distribution to points that grow and then shrink
# like the lightning bolts of the mighty zeus
# phase 1b -- shrink plot to the upper left corner, then make three copies into the other three corners, specifying each copy's nature before shrinking it (and having bodybuilders poke out)
# phase 2a -- fade in the different model specifications, as well as the prior distributions for each model's parameters in the appropriate subplots, zooming out as needed
# phase 2b -- pan a line across the plot, w/ either Bayes' head (+ the Stan logo) or Newton's head (+ the R logo) above (maybe give newton googly eyes?)
# ORRRR put arnold bayes and... 
# shifting prior of model parameters into posterior, and plot 89% HPDIs and 89% posterior predictive intervals in shaded region
# shade area to the left of the line with a light grey, and color open points according to the PSIS / LPPD score
# in total, have 4 models: 1) true model, flat priors 2) true model, shrinkage priors 3) cubic splines, shrinkage priors 4) true model, bfgs MLE w/ bootstrap proportions
# fade in gold star on winner, fade everything out

#### describe models in Stancode ####

### Stan code for model
logistic_model_sc <- "
data{
    int count[128];
    vector[128] time;
}
parameters{
    real logK_std;
    real<lower=0> r;
    real<lower=0> theta;
    real<lower=0, upper = 1> L0_prop;
}
transformed parameters{
    real<lower = 0> K = exp(log(10) + logK_std * 3 * log(10)); //non-centered parameterization
    real<lower = 0> L0 = K * L0_prop;;
}
model{
    vector[128] lambda;
    vector[128] alpha;
    vector[128] beta;
    theta ~ exponential( 1 );
    r ~ exponential( 2 );
    L0_prop ~ beta(1, 100);
    logK_std ~ std_normal();
    for ( i in 1:128 ) {
        lambda[i] = K/(1 + ((K - L0)/L0) * exp(-r * time[i]));
        alpha[i] = lambda[i] / theta;
        beta[i] = 1 / theta;
    }
    count ~ neg_binomial( alpha , beta );
}
generated quantities{

}
"

logistic_model_sc_flatpriors <- "
data{
    int count[128];
    vector[128] time;
}
parameters{
    real logK_std;
    real<lower=0> r_std;
    real<lower=0> theta_std;
    real<lower=0, upper = 1> L0_prop;
}
transformed parameters{
    real<lower = 0> K = exp(log(10) + logK_std * 50 * log(10)); //non-centered parameterization
    real<lower = 0> L0 = K * L0_prop;
    real<lower=0> r = r_std / 0.1; //non-centered parameterization
    real<lower=0> theta = theta_std / 0.01; //non-centered parameterization
}
model{
    vector[128] lambda;
    vector[128] alpha;
    vector[128] beta;
    theta_std ~ exponential( 1 );
    r_std ~ exponential( 1 );
    L0_prop ~ beta(1, 1);
    logK_std ~ std_normal();
    for ( i in 1:128 ) {
        lambda[i] = K/(1 + ((K - L0)/L0) * exp(-r * time[i]));
        alpha[i] = lambda[i] / theta;
        beta[i] = 1 / theta;
    }
    count ~ neg_binomial( alpha , beta );
}
generated quantities{

}
"


#### fit models ####

logistic_model_sc <- gsub(x = logistic_model_sc, pattern = "128", replacement = round(n_obs * prop_to_sample))
logistic_model_sc_flatpriors <- gsub(x = logistic_model_sc_flatpriors, pattern = "128", replacement = round(n_obs * prop_to_sample))
n_chains <- 4
initf <- function(chain_id = 1) {list(L0_prop = 0.1, r = 1, theta = 5, logK = 5, alpha = chain_id)}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

fit_regpriors <- stan(model_code = logistic_model_sc, data = d, chains = n_chains, warmup = 0.75E3, iter = 1E4, cores = 4, thin = 3, control = list(adapt_delta = 0.9))
fit_flatpriors <- stan(model_code = logistic_model_sc_flatpriors, data = d, chains = n_chains, warmup = 1.5E4, 
                       iter = 2E4, cores = 4, init = init_ll, thin = 3, control = list(adapt_delta = 0.95))

#### visualize fitted model results ####

for(fit in c(fit_regpriors, fit_flatpriors)){
  
logism <- extract.samples(fit)
thin <- 2
logism_expected_values <- sapply(1:(length(logism$logK) / thin), function(iter) 
  P(K = logism$K[iter*thin], r = logism$r[iter*thin], t = t, L0 = logism$L0[iter*thin]))
logism_89HPDI <- apply(logism_expected_values, 1, HPDI)

logism_sample_rates <- sapply(1:(length(logism$logK) / thin), function(iter) 
  rgamma(n = length(t), shape = logism_expected_values[,iter] / logism$theta[iter*thin], scale = logism$theta[iter*thin]))
logism_89HPDI_rates <- apply(logism_sample_rates, 1, HPDI)

logism_sample_obs <- sapply(1:(length(logism$logK) / thin), function(iter) 
  rpois(n = length(t), lambda = logism_sample_rates[,iter]))
logism_89HPDI_obs <- apply(logism_sample_obs, 1, HPDI)


plot(t, rates, type = "l", ylim = c(0,5000), lwd = 2, col = "darkred")
points(sampleTimes, sampleObs)

col_lambdaHPDI <- "#00416c"
col_sampHPDI <- "#c24e00"

lines(t, logism_89HPDI[1,], col = col_lambdaHPDI)
lines(t, logism_89HPDI[2,], col = col_lambdaHPDI)
polygon(y = c(logism_89HPDI[2,], rev(logism_89HPDI[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_lambdaHPDI, 0.5), border = NA)

# lines(t, logism_89HPDI_rates[1,], col = col_lambdaHPDI)
# lines(t, logism_89HPDI_rates[2,], col = col_lambdaHPDI)
# polygon(y = c(logism_89HPDI_rates[2,], rev(logism_89HPDI_rates[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_lambdaHPDI, 0.5))

lines(t, logism_89HPDI_obs[1,], col = col_sampHPDI)
lines(t, logism_89HPDI_obs[2,], col = col_sampHPDI)
polygon(y = c(logism_89HPDI_obs[2,], rev(logism_89HPDI[2,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)
polygon(y = c(logism_89HPDI[1,], rev(logism_89HPDI_obs[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)

}


#### fit all replicate models from 2 obs to 128 ####

#first the regularizing priors
do_regpriors <- F
n_chains <- 4
min_req_ess <- 2E3
n_acceptable_divergent_transitions <- 0
max_acceptable_rhat <- 1.001
if(do_regpriors){
  
  for(i in 5:2){
    
    fileout <- paste0("~/Documents/logistic_processes/just_growth/regpriors_", i, "_obs")
    
    #check if previous fit exists and load if so, to avoid recompiling
    if(file.exists(fileout)){
      load(fileout)
    }
    
    #get data subset
    dsub <- d
    dsub$time <- dsub$time[1:i]
    dsub$count <- dsub$count[1:i]
    
    #get model code
    logistic_model_sc_sub <- gsub(x = logistic_model_sc, pattern = n_obs, replacement = i)
    
    #specify fit requirements
    min_ess <- 0
    n_divergent_transitions <- 100
    max_rhat <- 10
    n_warmup <- 0.75E3
    n_iter <- 1E4
    adapt_delta = 0.8
    thin = 2
    while(min_ess < min_req_ess | n_divergent_transitions > n_acceptable_divergent_transitions | max_rhat > max_acceptable_rhat){
      
      cat(paste0(i, " "))
      
      #fit model
      if(!exists("fit")){
        fit <- stan(model_code = logistic_model_sc_sub, data = dsub, chains = n_chains, warmup = n_warmup, 
                    iter = n_iter, cores = 4, thin = thin, control = list(adapt_delta = adapt_delta))
      } else{
        fit <- stan(fit = fit, data = dsub, chains = n_chains, warmup = n_warmup, 
                    iter = n_iter, cores = 4, thin = thin, control = list(adapt_delta = adapt_delta))
      }
      
      #do automated diagnostics
      sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
      n_divergent_transitions <- max(sapply(sampler_params, function(x) max(x[, "divergent__"])))
      max_rhat <- max(summary(fit)$summary[, "Rhat"])
      min_ess <- min(summary(fit)$summary[, "n_eff"])
      
      #make necessary adjustments
      if(min_ess < min_req_ess | max_rhat > max_acceptable_rhat){
        n_warmup <- n_warmup * 2
        n_iter <- n_iter * 2
        thin <- thin * 2
        print(paste0("ess = ", min_ess))
        print(paste0("rhat = ", max_rhat))
      } 
      if(n_divergent_transitions > n_acceptable_divergent_transitions){
        adapt_delta <- mean(c(adapt_delta, 1))
        print(paste0("# divergent transitions: ", n_divergent_transitions))
      } 
      
    }
    
    #save stanfit object
    save(fit, file = fileout)
    rm(fit)
    
  }
  
}

#now the flat priors
do_flatpriors <- T
n_chains <- 4
min_req_ess <- 2E3
n_acceptable_divergent_transitions <- 1
max_acceptable_rhat <- 1.001
initf <- function(chain_id = 1) {list(L0_prop = 0.05, r = 1, theta = 5, logK = 5, alpha = chain_id)}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))
if(do_flatpriors){
  
  for(i in 22:2){
    
    fileout <- paste0("~/Documents/logistic_processes/just_growth/flatpriors_", i, "_obs")
    
    #check if previous fit exists and load if so, to avoid recompiling
    if(file.exists(fileout)){
      load(fileout)
    }
    
    #get data subset
    dsub <- d
    dsub$time <- dsub$time[1:i]
    dsub$count <- dsub$count[1:i]
    
    #get model code
    logistic_model_sc_flatpriors_sub <- gsub(x = logistic_model_sc_flatpriors, pattern = n_obs, replacement = i)
    
    #specify fit requirements
    min_ess <- 0
    n_divergent_transitions <- 100
    max_rhat <- 10
    n_warmup <- 1.5E3
    n_iter <- 2E4
    adapt_delta = 0.9
    thin = 2
    while(min_ess < min_req_ess | n_divergent_transitions > n_acceptable_divergent_transitions | max_rhat > max_acceptable_rhat){
      
      cat(paste0(i, " "))
      
      #fit model
      if(!exists("fit")){
        fit <- stan(model_code = logistic_model_sc_flatpriors_sub, data = dsub, chains = n_chains, warmup = n_warmup, 
                    iter = n_iter, cores = 4, thin = thin, init = init_ll, control = list(adapt_delta = adapt_delta))
      } else{
        fit <- stan(fit = fit, data = dsub, chains = n_chains, warmup = n_warmup, init = init_ll,
                    iter = n_iter, cores = 4, thin = thin, control = list(adapt_delta = adapt_delta))
      }
      
      #do automated diagnostics
      sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
      n_divergent_transitions <- max(sapply(sampler_params, function(x) max(x[, "divergent__"])))
      max_rhat <- max(summary(fit)$summary[, "Rhat"])
      min_ess <- min(summary(fit)$summary[, "n_eff"])
      
      #make necessary adjustments
      if(min_ess < min_req_ess | max_rhat > max_acceptable_rhat){
        n_warmup <- n_warmup * 2
        n_iter <- n_iter * 2
        thin <- thin * 2
        print(paste0("ess = ", min_ess))
        print(paste0("rhat = ", max_rhat))
      } 
      if(n_divergent_transitions > n_acceptable_divergent_transitions){
        adapt_delta <- mean(c(adapt_delta, 1))
        print(paste0("# divergent transitions: ", n_divergent_transitions))
      } 
      
    }
    
    #save stanfit object
    save(fit, file = fileout)
    rm(fit)
    
  }
  
}

#### now use optim to find the MLE #####

data <- as.data.frame(d)
#par is c(K, L0_prop, r, theta)
gampois_logisticgrowth <- function(data, par) {

    #undo params
    K <- par[1]
    L0_prop <- par[2]
    r <- par[3]
    theta <- par[4]

    #impose bounds on Nelder-Mead muahaha
    if(L0_prop < 0 | L0_prop > 1){return(Inf)}

    #model
    L0 = K * L0_prop
    lambda = K /(1 + ((K - L0)/L0) * exp(-r * data$time))
    alpha = lambda / theta
    beta = 1 / theta
    nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
    # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
    # print(c(nloglik, nloglik2)) #these evaluate to the same value
    
    if(is.na(nloglik)){
      cat(paste0("\nlog-likelihood is ", nloglik, "\n"))
      cat(paste0("\nr is ", r, "\n"))
      cat(paste0("\ntheta is ", theta, "\n"))
      cat(paste0("\nK is ", K, "\n"))
      cat(paste0("\nL0 is ", L0, "\n"))
    }

    #penalized likelihood?
    nloglik <- nloglik + -dnorm(log10(K), mean = 1, sd = 3, log = T)
    nloglik <- nloglik + -dexp(theta, 1, log = T)
    nloglik <- nloglik + -dexp(r, 2, log = T)
    nloglik <- nloglik + -dbeta(x = L0_prop, shape1 = 1, shape2 = 100, log = T)

    return(nloglik)

}
gampois_logisticgrowth(data, par = c(100, 0.1, 1, 1))
est_all <- optimx::optimx(par = c(10, 0.5, 1, 1), fn = gampois_logisticgrowth, data = data, 
             method = "nlminb", lower = c(1E-6, 1E-6, 1E-6, 1E-6), upper = c(Inf, 1-1E-6, Inf, Inf), control = list(maxit = 1E4, trace = 2))
# est <- optim(par = c(1000, 0.01, 1, 1), fn = gampois_logisticgrowth, data = data,
#              method = "Nelder-Mead", control = list(trace = 1))


est <- as.list(as.numeric(est_all[names(est_all)[grep(names(est_all), pattern = "p")]]))
names(est) <- c("K", "L0_prop", "r", "theta")
est$L0 <- est$K * est$L0_prop
est_L <- P(K = est$K, r = est$r, t = t, L0 = est$L0)

plot(t, rates, type = "l", ylim = c(0,5000), lwd = 2, col = "darkred")
points(sampleTimes, sampleObs)
lines(t, est_L, col = "#009ECE", lwd = 2)
lines(t, rates, col = "red", lwd = 2)
legend(x = "topleft", col=c("red", "#009ECE"), legend = c("true model", "estimated model"), lwd = 2)


#### show samples from the prior ####

#regularizing prior
nsamp <- 10000
Ks <- 10^rnorm(nsamp, mean = 1, sd = 3)
thetas <- rexp(nsamp, 1)
rs <- rexp(nsamp, 2)
L0_props <- rbeta(nsamp, shape1 = 1, shape2 = 100)
L0s <- Ks * L0_props
ratess <- t(sapply(1:nsamp, function(x) P(K = Ks[x], r = rs[x], t = t, L0 = L0s[x])))
prior_89HPDI <- apply(ratess, 2, HPDI)
prior_sample_rates <- sapply(1:nsamp, function(iter) 
  rgamma(n = length(t), shape = ratess[iter,] / thetas[iter], scale = thetas[iter]))
prior_sample_obs <- sapply(1:nsamp, function(iter) rpois(n = length(t), lambda = prior_sample_rates[,iter]))
prior_89HPDI_obs <- apply(prior_sample_obs, 1, HPDI)

plot(t, rates, type = "l", ylim = c(0,max(prior_89HPDI_obs)), lwd = 2, col = "darkred")
points(sampleTimes, sampleObs)
nlines_toplot <- 4E3
for(i in 1:nlines_toplot){lines(t, ratess[i * trunc(nsamp/nlines_toplot),], col = grDevices::adjustcolor("black", 0.35), lwd = 1)}

col_lambdaHPDI <- "#00416c"
col_sampHPDI <- "#c24e00"
lines(t, prior_89HPDI[1,], col = col_lambdaHPDI)
lines(t, prior_89HPDI[2,], col = col_lambdaHPDI)
polygon(y = c(prior_89HPDI[2,], rev(prior_89HPDI[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_lambdaHPDI, 0.5), border = NA)
lines(t, prior_89HPDI_obs[1,], col = col_sampHPDI)
lines(t, prior_89HPDI_obs[2,], col = col_sampHPDI)
polygon(y = c(prior_89HPDI_obs[2,], rev(prior_89HPDI[2,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)
polygon(y = c(prior_89HPDI[1,], rev(prior_89HPDI_obs[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)

#flat prior
nsamp <- 1E6
Ks <- 10^rnorm(nsamp, mean = 1, sd = 50)
thetas <- rexp(nsamp, 0.01)
rs <- rexp(nsamp, 0.1)
L0_props <- rbeta(nsamp, shape1 = 1, shape2 = 1)
L0s <- Ks * L0_props
ratess <- t(sapply(1:nsamp, function(x) P(K = Ks[x], r = rs[x], t = t, L0 = L0s[x])))
prior_89HPDI <- apply(ratess, 2, HPDI)
prior_sample_rates <- sapply(1:nsamp, function(iter) 
  rgamma(n = length(t), shape = ratess[iter,] / thetas[iter], scale = thetas[iter]))
prior_sample_obs <- sapply(1:nsamp, function(iter) rpois(n = length(t), lambda = prior_sample_rates[,iter]))
prior_89HPDI_obs <- apply(prior_sample_obs, 1, HPDI)

plot(t, rates, type = "l", ylim = c(0,max(prior_89HPDI_obs)), lwd = 2, col = "darkred")
points(sampleTimes, sampleObs)
nlines_toplot <- 4E3
for(i in 1:nlines_toplot){lines(t, ratess[i * trunc(nsamp/nlines_toplot),], col = grDevices::adjustcolor("black", 0.35), lwd = 1)}

col_lambdaHPDI <- "#00416c"
col_sampHPDI <- "#c24e00"
lines(t, prior_89HPDI[1,], col = col_lambdaHPDI)
lines(t, prior_89HPDI[2,], col = col_lambdaHPDI)
polygon(y = c(prior_89HPDI[2,], rev(prior_89HPDI[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_lambdaHPDI, 0.5), border = NA)
lines(t, prior_89HPDI_obs[1,], col = col_sampHPDI)
lines(t, prior_89HPDI_obs[2,], col = col_sampHPDI)
polygon(y = c(prior_89HPDI_obs[2,], rev(prior_89HPDI[2,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)
polygon(y = c(prior_89HPDI[1,], rev(prior_89HPDI_obs[1,])), x = c(t, rev(t)), col = grDevices::adjustcolor(col_sampHPDI, 0.5), border = NA)

#### do some bootstrapping & nonlinear optimization ####

data <- as.data.frame(d)
#par is c(K, L0_prop, r, theta)
gampois_logisticgrowth_penlik <- function(data, par) {
  
  #undo params
  K <- par[1]
  L0_prop <- par[2]
  r <- par[3]
  theta <- par[4]
  
  #impose bounds on unbounded algos muahaha
  if(L0_prop < 0 | L0_prop > 1 | K < 0 | theta < 0){return(Inf)}
  
  #model
  L0 = K * L0_prop
  lambda = K /(1 + ((K - L0)/L0) * exp(-r * data$time))
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood is ", nloglik, "\n"))
    cat(paste0("\nr is ", r, "\n"))
    cat(paste0("\ntheta is ", theta, "\n"))
    cat(paste0("\nK is ", K, "\n"))
    cat(paste0("\nL0 is ", L0, "\n"))
  }
  
  #penalized likelihood?
  nloglik <- nloglik + -dnorm(log10(K), mean = 1, sd = 3, log = T)
  nloglik <- nloglik + -dexp(theta, 1, log = T)
  nloglik <- nloglik + -dexp(r, 2, log = T)
  nloglik <- nloglik + -dbeta(x = L0_prop, shape1 = 1, shape2 = 100, log = T)
  
  return(nloglik)
  
}
gampois_logisticgrowth_nopenlik <- function(data, par) {
  
  #undo params
  K <- par[1]
  L0_prop <- par[2]
  r <- par[3]
  theta <- par[4]
  
  #impose bounds on unbounded algos muahaha
  if(L0_prop < 0 | L0_prop > 1 | K < 0 | theta < 0){return(Inf)}
  
  #model
  L0 = K * L0_prop
  lambda = K /(1 + ((K - L0)/L0) * exp(-r * data$time))
  alpha = lambda / theta
  beta = 1 / theta
  nloglik <-  -sum(extraDistr::dgpois(data$count, shape = alpha, rate = beta, log = T))
  # nloglik2 <-  -sum(sads::dpoig(data$count, shape = alpha, rate = beta, log = T, frac = 1))
  # print(c(nloglik, nloglik2)) #these evaluate to the same value
  
  if(is.na(nloglik)){
    cat(paste0("\nlog-likelihood is ", nloglik, "\n"))
    cat(paste0("\nr is ", r, "\n"))
    cat(paste0("\ntheta is ", theta, "\n"))
    cat(paste0("\nK is ", K, "\n"))
    cat(paste0("\nL0 is ", L0, "\n"))
  }
  
  return(nloglik)
  
}

library(foreach)
library(parallel)
library(doParallel)
do_bootstrap <- T
n_bootstrap_replicates <- 2E3
if(do_bootstrap){
  
  if(!exists("cl")){
    cl <- makeCluster(8, outfile="")
    registerDoParallel(cl)
  }
  
  getDoParWorkers()
  
  foreach(i=128:1, .packages = c("optimx")) %dopar% {
    
    for(lik in c("nopenlik", "penlik")){
    
      fileout <- paste0("~/Documents/logistic_processes/just_growth/mle_bootstrap_", lik, "_", i, "_obs")
      
      #get data subset
      dsub <- d
      dsub$time <- dsub$time[1:i]
      dsub$count <- dsub$count[1:i]
      data <- as.data.frame(dsub)
      
      #resample data w/ replacement
      data_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) sample(x = 1:i, size = i, replace = T))
      data_bootstrap <- lapply(1:n_bootstrap_replicates, function(bsr) data[data_bootstrap[[bsr]],])
      
      #initialize output object
      fit <- list(MLE = NA, BSREPS = matrix(NA, nrow = n_bootstrap_replicates, ncol= 5))
      
      #find MLE
      est_all <- optimx::optimx(par = c(10, 0.5, 1, 1), 
                                fn = ifelse(lik == "penlik", gampois_logisticgrowth_penlik, gampois_logisticgrowth_nopenlik), 
                                data = data, method = "nlm", lower = c(1E-6, 1E-6, 1E-6, 1E-6), upper = c(Inf, 1-1E-6, Inf, Inf), 
                                control = list(maxit = 1E4, trace = 0))
      est <- as.list(as.numeric(est_all[names(est_all)[grep(names(est_all), pattern = "p")]]))
      names(est) <- c("K", "L0_prop", "r", "theta")
      est$L0 <- est$K * est$L0_prop
      fit$MLE <- unlist(est)
      
      #find MLE for all bootstrap replicates
      for(j in 1:n_bootstrap_replicates){
        cat(paste0("(", i, ", ", j, ") " ))
        est_all <- optimx::optimx(par = c(10, 0.5, 1, 1), fn = ifelse(lik == "penlik", gampois_logisticgrowth_penlik, gampois_logisticgrowth_nopenlik), 
                                  data = data_bootstrap[[j]], method = "nlm", lower = c(1E-6, 1E-6, 1E-6, 1E-6), upper = c(Inf, 1-1E-6, Inf, Inf), 
                                  control = list(maxit = 1E4, trace = 0))
        est <- as.list(as.numeric(est_all[names(est_all)[grep(names(est_all), pattern = "p")]]))
        names(est) <- c("K", "L0_prop", "r", "theta")
        est$L0 <- est$K * est$L0_prop
        fit$BSREPS[j,] <- unlist(est)
      }
      
      #save bootstrap object
      save(fit, file = fileout)
      rm(fit)
    
    }
  
  
  }

  stopCluster(cl)
  rm(cl)
  
}