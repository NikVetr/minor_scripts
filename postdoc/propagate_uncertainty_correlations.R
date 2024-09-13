library(cmdstanr)
library(posterior)
par_base <- par(no.readonly = TRUE)

#set a seed for reproducibility
set.seed <- 123

#specify model with uncertainty
model_index <- 1
model_loc <- list("~/scripts/minor_scripts/postdoc/correlation_uncertainty.stan", #1
                  "~/scripts/minor_scripts/postdoc/correlation_uncertainty_horseshoe-unpooled.stan", #2
                  "~/scripts/minor_scripts/postdoc/correlation_uncertainty_horseshoe-partially-pooled.stan", #3
                  "~/scripts/minor_scripts/postdoc/correlation_uncertainty_horseshoe-pooled.stan", #4
                  "~/scripts/minor_scripts/postdoc/correlation_uncertainty_logit-rescale.stan", #5
                  "~/scripts/minor_scripts/postdoc/correlation_uncertainty_bivariate-laplace.stan", #6
                  "", #7
                  "", #8
                  "", #9
                  "", #10
                  "", #11
                  "", #12
                  "", #13
                  "" #14
)[[model_index]]
model_string <-  paste0(readLines(model_loc), collapse = "\n")

if(grepl("bivariate-laplace", model_loc)){
    bivariate_laplace = T;
} else {
    bivariate_laplace = F;    
}


#can maybe add a further term to the model to represent distribution of estimated errors?


#### simulate input ####

#simulate coefficients
r <- 0.75 #true correlation between coefficients
R <- diag(2) * (1-r) + r #corresponding correlation matrix
p <- 2E3 #total number of samples
x <- matrix(rnorm(p*2), ncol = 2) %*% chol(R) #sample true coefficients

#modify data according to model

if(grepl("horseshoe", model_loc)){
    prop_null <- 0
    if(prop_null != 0){
        # null_inds <- sample(1:(2*n), ceiling(n*prop_null*2))
        # x[null_inds] <- x[null_inds] / 10
        
        n_null <- ceiling(p*prop_null)
        null_inds <- sample(1:p, n_null)
        which_null <- 1 + rbinom(n_null, 2, prob = 0.5)
        x[null_inds[which_null == 2],] <- x[null_inds[which_null == 2],] / 10
        x[null_inds[which_null == 1], 1] <- x[null_inds[which_null == 1], 1] / 10
        x[null_inds[which_null == 3], 2] <- x[null_inds[which_null == 3], 2] / 10
    }
}

if(bivariate_laplace){
    exp_rv <- rexp(p, 1)
    x[,1] <- x[,1] * exp_rv
    x[,2] <- x[,2] * exp_rv
}

#simulate data and fit OLS model
n <- c(20, 5) #define sample size across dimensions
e_sd_sd <- c(1, 1) #define differential power across two dimensions? or just get from sample size
e_sd <- matrix(rexp(p*2), ncol = 2) %*% diag(e_sd_sd) #sample element-wise error
sim_and_fit_lm <- function(b, err_sd, nobs){
    asim <- rnorm(1)
    xsim <- rnorm(nobs)
    esim <- rnorm(n = nobs, mean = 0, sd = err_sd)
    ysim <- asim + xsim * b + esim
    fit <- lm(ysim ~ xsim)
    summary(fit)$coefficients[2,1:2]
}

fits_1 <- do.call(rbind, lapply(1:p, function(i){
    sim_and_fit_lm(b = x[i,1], err_sd = e_sd[i,1], nobs = n[1])
}))
fits_2 <- do.call(rbind, lapply(1:p, function(i){
    sim_and_fit_lm(b = x[i,2], err_sd = e_sd[i,2], nobs = n[2])
}))

x_err <- cbind(fits_1[,"Estimate"], fits_2[,"Estimate"])
sd_x_err <- cbind(fits_1[,"Std. Error"], fits_2[,"Std. Error"])

#plot these
par(mfrow = c(1,2), mar = c(5,5,5,1))
plot(x, xlab = "true coefficient 1", 
     ylab = "true coefficient 2", 
     main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                  round(cor(x)[1,2], 3))),
     col = adjustcolor(1, 0.3), pch = 19)
abline(0,1, col = 2, lty = 2, lwd = 2)

plot(x_err, xlab = "OLS estimated coefficient 1", 
     ylab = "OLS estimated coefficient 2", 
     main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                  round(cor(x_err)[1,2], 3))),
     col = adjustcolor(1, 0.3), pch = 19)
abline(0,1, col = 2, lty = 2, lwd = 2)

#evaluate sample correlations
cor(x)[1,2] #should be approximately r
cor(x_err)[1,2] #should be regressed from r towards 0


#### preprocess model ####

#pack data into a list
dat <- list(p=p, x_err=x_err, sd_x_err=sd_x_err)

#or for the horseshoe version
if(grepl("horseshoe", model_loc)){
    dat <- list(p=p, x_err=x_err, sd_x_err=sd_x_err, n = mean(n),
                p0 = c(rep(p, ifelse(grepl("unpooled", model_loc), 2, 1))) * (1-prop_null))
}

#or for the logit rescale version
if(grepl("logit-rescale", model_loc)){
    dat <- list(p=p, x_err=x_err, sd_x_err=sd_x_err, 
                p0 = c(rep(p, ifelse(grepl("unpooled", model_loc), 2, 1))) * (1-prop_null))
}

#or for the bivariate laplace
if(bivariate_laplace){
    dat <- list(n=p, x_err=x_err, sd_x_err=sd_x_err)    
}

#### fit model ####
mod <- cmdstan_model(model_loc)
fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                  data = dat, adapt_delta = 0.9, parallel_chains = 4,
                  refresh = 10, max_treedepth = 10, 
                  thin = 1, init = 0.1)

#### mcmc diagnostics ####

#check convergence
summ <- fit$summary()
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])

#extract samples and inspect
samps <- data.frame(as_draws_df(fit$draws()))
samps_r <- data.frame(as_draws_df(fit$draws("r")))
samps_x <- data.frame(as_draws_df(fit$draws("x")))

#### plotting ####
cols <- adjustcolor(c(2,4), 0.5)
cols_vec <- rep(cols, each = p)
shuffi <- sample(1:(p*2))

#plot marginal posterior
range_breaks <- range(c(samps_r$r, r, cor(x_err)[1,2]))
breaks <- floor(range_breaks[1]*50):ceiling(range_breaks[2]*50)/50
hist(samps_r$r, breaks = breaks, 
     main = "posterior distribution\nof correlation", freq = F, xlab = "correlation")

#label lines
abline(v = r, col = "red", lwd = 3)
text(x = r, y = par("usr")[4], pos = 3, labels = "true\ncorr.", xpd = NA, col = "red")
abline(v = cor(x_err)[1,2], col = "blue", lwd = 3)
text(x = cor(x_err)[1,2], y = par("usr")[4], 
     pos = 3, labels = "sample\ncorr.", xpd = NA, col = "blue")

#look at estimates of x
est_x_pmeans <- apply(samps_x, 2, mean)
x_pmeans <- cbind(est_x_pmeans[1:p], est_x_pmeans[(1+p):(2*p)])
est_x_psds <- apply(samps_x, 2, sd)
x_psds <- cbind(est_x_psds[1:p], est_x_psds[(1+p):(2*p)])

#plot means
par(mfrow = c(1,2), mar = c(5,5,5,1))
plot(x[shuffi], x_err[shuffi], xlab = "true coefficient", 
     ylab = latex2exp::TeX("\\textbf{OLS estimate} of coefficient"), 
     main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                             round(cor(cbind(c(x), c(x_err)))[1,2], 3))),
     col = cols_vec[shuffi], pch = 19)
abline(0,1, col = 1, lty = 2, lwd = 2)

plot(x[shuffi], x_pmeans[shuffi], xlab = "true coefficient", 
     ylab = latex2exp::TeX("\\textbf{posterior mean} of coefficient"),
     main = latex2exp::TeX(paste0("Pearson's \\textit{r} ≈ ", 
                                               round(cor(cbind(c(x), c(x_pmeans)))[1,2], 3))),
     col = cols_vec[shuffi], pch = 19)
abline(0,1, col = 1, lty = 2, lwd = 2)

#plot sds
par(par_base)

#scatterplot
plot(log10(sd_x_err[shuffi]), log10(x_psds[shuffi]), pch = 19, 
     col = cols_vec[shuffi],
     xlab = latex2exp::TeX("log$_{10}$(OLS Standard Error)"), 
     ylab = latex2exp::TeX("log$_{10}$Posterior SD)"))
legend("topleft", legend=c("Var 1", "Var 2"), 
    col = cols, pch = 19)
abline(0,1,lwd=2,lty=2)

#get histogram data
sd_diffs <- log10(sd_x_err) - log10(x_psds) #positive means posterior sd is smaller
breaks <- seq(min(sd_diffs), max(sd_diffs), length.out = 20)
histall <- hist(sd_diffs, breaks = breaks, plot=FALSE)
hist1 <- hist(sd_diffs[,1], breaks = breaks, plot=FALSE)
hist2 <- hist(sd_diffs[,2], breaks = breaks, plot=FALSE)

#munge to pseudolog for frequency
hist1$counts <- log10(hist1$counts + 1)
hist2$counts <- log10(hist2$counts + 1)
histall$counts <- log10(histall$counts + 1)

# Plot the stacked histogram
par(fig=c(0.5, 0.95, 0.25, 0.5), mar = c(0,0,0,0), new=TRUE)
plot(hist1, col = cols[1], xlab="", ylab="", main = "", xpd = NA, cex.axis=0.7)
mtext(1, text = "Difference", line = 2, cex = 0.7)
mtext(2, text = latex2exp::TeX("log$_{10}$(Frequency)"), line = 2, cex = 0.7)
plot(hist2, col=cols[2], add=TRUE)

#reset par
par(par_base)

