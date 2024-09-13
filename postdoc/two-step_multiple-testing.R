#simulation function for first step of the analysis
sim_analysis_step1 <- function(b_x1, b_x2, n = 1E2, sigmas = c(1,1,1)){
  
  #sim inputs
  x1 <- rnorm(n) * sigmas[1]
  x2 <- rnorm(n) * sigmas[2]
  d <- data.frame(x = c(x1, x2), f = rep(1:2, each = n))
  
  #sim outputs
  e <- rnorm(n*2) * sigmas[3]
  d$y <- d$x * c(b_x1, b_x2)[d$f] + e

  #run first step of 2 step regression
  fit_2step_A <- lm(d$y[d$f==1] ~ 1 + x1)
  fit_2step_B <- lm(d$y[d$f==2] ~ 1 + x2)
  pvals_2step <- c(A = summary(fit_2step_A)$coefficients["x1", "Pr(>|t|)"],
                   B = summary(fit_2step_B)$coefficients["x2", "Pr(>|t|)"])
  
  #return pvals and data
  return(list(pvals = pvals_2step,
              d = d))
  
}

#simulate first stage
nsim <- 1E4
b <- 0.1
fits <- parallel::mclapply(1:nsim, function(i) sim_analysis_step1(b_x1 = b, b_x2 = b), mc.cores = 4)
pvals_step1 <- data.frame(do.call(rbind, lapply(fits, function(x) x[["pvals"]])))
q_thresh <- 0.05
qvals_step1 <- data.frame(A = p.adjust(pvals_step1$A, method = "BH"),
                          B = p.adjust(pvals_step1$B, method = "BH"))
sig_tests_step1 <- apply(qvals_step1 < q_thresh, 1, any)
sig_fits_step1 <- fits[sig_tests_step1]

#fit interaction model to filtered results from first stage
pvals_step2 <- unlist(parallel::mclapply(sig_fits_step1, function(x){
  d <- x[["d"]]
  fit <- lm(y ~ 1 + x + f + x*f, d)
  return(summary(fit)$coefficients["x:f", "Pr(>|t|)"])
}, mc.cores = 4))

#plot results
par(mfrow = c(2,1))
hist(pvals_step2, breaks = 0:40/40, 
     main = "nominal interaction pvals (true interaction effect = 0)", 
     xlab = "p-value")
hist(p.adjust(pvals_step2, method = "BH"), breaks = 0:40/40, 
     main = "fdr-adjusted interaction qvals (true interaction effect = 0)", 
     xlab = "p-value")

#alternatively, filter only for main effects in combined analysis
sim_analysis_step1_combined <- function(b_x1, b_x2, n = 1E2, sigmas = c(1,1,1)){
  
  #sim inputs
  x1 <- rnorm(n) * sigmas[1]
  x2 <- rnorm(n) * sigmas[2]
  d <- data.frame(x = c(x1, x2), f = rep(1:2, each = n))
  
  #sim outputs
  e <- rnorm(n*2) * sigmas[3]
  d$y <- d$x * c(b_x1, b_x2)[d$f] + e
  
  #run first step of 2 step regression
  fit_2step_comb <- lm(d$y ~ 1 + d$x)
  
  #return pvals and data
  return(list(pvals = summary(fit_2step_comb)$coefficients["d$x", "Pr(>|t|)"],
              d = d))
  
}

#simulate first stage
fits_comb <- parallel::mclapply(1:nsim, function(i) sim_analysis_step1_combined(b_x1 = b, b_x2 = b), mc.cores = 4)
pvals_step1_comb <- unlist(lapply(fits_comb, function(x) x[["pvals"]]))
qvals_step1_comb <- p.adjust(pvals_step1_comb, method = "BH")
sig_tests_step1_comb <- qvals_step1_comb < q_thresh
sig_fits_step1_comb <- fits_comb[sig_tests_step1_comb]

#fit interaction model to filtered results from first stage
pvals_step2_comb <- unlist(parallel::mclapply(sig_fits_step1_comb, function(x){
  d <- x[["d"]]
  fit <- lm(y ~ 1 + x + f + x*f, d)
  return(summary(fit)$coefficients["x:f", "Pr(>|t|)"])
}, mc.cores = 4))

#plot results
par(mfrow = c(2,1))
hist(pvals_step2_comb, breaks = 0:40/40, 
     main = "nominal interaction pvals (true interaction effect = 0)", 
     xlab = "p-value")
hist(p.adjust(pvals_step2_comb, method = "BH"), breaks = 0:40/40, 
     main = "fdr-adjusted interaction qvals (true interaction effect = 0)", 
     xlab = "p-value")
