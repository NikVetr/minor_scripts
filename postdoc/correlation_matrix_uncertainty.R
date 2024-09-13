library(cmdstanr)
library(posterior)
library(dplyr)
library(data.table)
source("~/repos/Stan2R/R/functions.R")

#specify a few functions
rlkj <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}

#### simulate data  ####

#generate high level params
rop <- 0.8
p <- 30
rs <- c(1, rop^(1:p))
R <- outer(1:p, 1:p, FUN = function(i, j, rs) rs[abs(i - j) + 1], rs = rs)

#or sample from lkj
R <- rlkj(p)

#get true eigendecomp
true_eigen <- eigen(R)

#sample data
n <- 10
x <- matrix(rnorm(n*p), n, p)
if(n>p){
  # x <- x %*% solve(chol(cov(x))) %*% chol(R)
  x <- x %*% chol(R)
} else {
  x <- x %*% chol(R) #whitening will just put 0s on the diag 
}

#evaluate correlations quickly
s_R <- cor(x)
Rpd <- as.matrix(Matrix::nearPD(s_R, corr = T)$mat)
ev_Rpd <- eigen(Rpd)$values

#evaluate pairwise correlations
combos <- t(apply(expand.grid(1:p, 1:p), 1, sort))
combos <- combos[!duplicated(combos),]
combos <- combos[!apply(combos, 1, function(x) x[1] == x[2]),]
pairwise_corrs <- data.frame(do.call(rbind, lapply(1:choose(p, 2), function(i){
  r <- cor(x[,combos[i,1]], x[,combos[i,2]])
  n <- nrow(x)
  z <- 0.5 * log((1 + r) / (1 - r))
  sigma_z <- 1 / sqrt(n - 3)
  return(c(z = z, sigma_z = sigma_z))
})))


#specify model with uncertainty
# model_loc <- "~/scripts/minor_scripts/postdoc/correlation_uncertainty.stan"
model_loc <- "~/scripts/minor_scripts/postdoc/correlation-matrix_from-bivariate-Fisher-z.stan"
model_string <-  paste0(readLines(model_loc), collapse = "\n")

#pack data into a list
dat <- list(p = p,
            np = choose(p,2),
            ip = combos,
            z = pairwise_corrs$z,
            sigma_z = pairwise_corrs$sigma_z)
#or for the horseshoe version

#### fit model ####
mod <- cmdstan_model(model_loc)
fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3,
                  data = dat, adapt_delta = 0.85, parallel_chains = 4,
                  refresh = 100, max_treedepth = 10, 
                  thin = 1, init = 0.01)

#check convergence
# summ <- fit$summary()
# print(summ[order(summ$ess_bulk),])
# print(summ[order(summ$rhat, decreasing = T),])

#extract samples and inspect
samps <- data.frame(as_draws_df(fit$draws()))
R_samps <- subset_samps("R", as.data.table(samps))
R_pmean <- apply(R_samps, 2, mean)
R_psd <- apply(R_samps, 2, sd)
Rpm <- munge_samps("R", rbind(R_pmean,R_pmean))[[1]]
pairs(cbind(R = R[upper.tri(R)], Rpm = Rpm[upper.tri(R)], Rpd = Rpd[upper.tri(R)]))
cor(Rpm[upper.tri(R)], R[upper.tri(R)])
cor(Rpd[upper.tri(R)], R[upper.tri(R)])

cor(Rpm[upper.tri(R)], R[upper.tri(R)])
cor(Rpd[upper.tri(R)], R[upper.tri(R)])

#### look at distribution of determinants, eigenvalues, etc ####

#for nearPD matrix
det_Rpd <- det(Rpd)
ev_Rpd <- eigen(Rpd)$values
evec_Rpd <- eigen(Rpd)$vectors

#for posterior sample
R_samps <- munge_samps("R", subset_samps("R", as.data.table(samps)))
det_samps <- sapply(R_samps, det)
eval_samps <- do.call(rbind, lapply(R_samps, function(foo) eigen(foo)$values))

find_angle <- function(vec1, vec2){
  rads <- acos(sum(vec1 * vec2))
  rads <- min(abs(rads), abs(pi - rads))
  return(rads / 2 / pi * 360)
}

evec_angle_samps <- do.call(rbind, lapply(R_samps, function(si){
  evs <- eigen(si)$vectors
  sapply(1:p, function(evi){
    find_angle(evs[,evi], true_eigen$vectors[,evi])
  })
}))

Rpd_angles <- sapply(1:p, function(evi){
  find_angle(evec_Rpd[,evi], true_eigen$vectors[,evi])
})

#plot evecs
#lower proportion is better
relative_goodness <- sapply(1:p, function(evi){ 
  mean(evec_angle_samps[,evi] > Rpd_angles[evi]) * 100
})
hist(relative_goodness, breaks = 0:10*10)
relative_goodness[1]

#### plotting EVs ####
hist(log10(det_samps), breaks = 100)
abline(v = log10(det(R)), lwd = 2, col = 3)
log10(det_Rpd)

bad_evs <- which(ev_samps <= 0, arr.ind = T)[,1]
ev_subsamps <- ev_samps
if(length(bad_evs) != 0){
  ev_subsamps <- ev_samps[-bad_evs,]
}
ev_quantiles <- apply(log10(ev_subsamps), 2, quantile, prob = c(0.005,0.05, 0.25, 0.75, 0.95, 0.995))
par(mar = c(6,6,2,2))
plot(log10(ev_Rpd), type = "l", lwd = 2, col = 2, ylim = range(c(ev_quantiles, log10(ev_Rpd))), 
     xlab = "Index of Eigenvalue", ylab = latex2exp::TeX("log$_{10}$(eigenvalue)"))
# lines(ev_quantiles[1,], lty = 3, col = adjustcolor(1, 0.5))
# lines(ev_quantiles[4,], lty = 3, col = adjustcolor(1, 0.5))
# lines(ev_quantiles[2,], lty = 3, col = adjustcolor(1, 0.5))
# lines(ev_quantiles[3,], lty = 3, col = adjustcolor(1, 0.5))
polygon(x = c(1:p,p:1,1), y = c(ev_quantiles["0.5%",], rev(ev_quantiles["99.5%",]), ev_quantiles["0.5%",1]), col = adjustcolor(1, 0.25), border = NA)
polygon(x = c(1:p,p:1,1), y = c(ev_quantiles["5%",], rev(ev_quantiles["95%",]), ev_quantiles["5%",1]), col = adjustcolor(1, 0.25), border = NA)
polygon(x = c(1:p,p:1,1), y = c(ev_quantiles["25%",], rev(ev_quantiles["75%",]), ev_quantiles["25%",1]), col = adjustcolor(1, 0.25), border = NA)
lines(log10(true_eigen$values), col = 5, type = "l", lwd = 3)
lines(log10(ev_Rpd), col = 2, type = "l", lwd = 3)
legend(x = par("usr")[1] + diff(par("usr")[1:2])/20, y = par("usr")[3] + diff(par("usr")[3:4])/2.5,
       legend = c("sample eigenvalues", 
                  "true eigenvalues", 
                  "eigenvalue CIs\n(99%, 90%, 50%)"),
       lwd = c(3,3,NA),
       col = c(2, 5, adjustcolor(1,0.4)),
       pch = c(NA, NA, 15),
       cex = 0.9, box.lwd = 0, pt.cex = 3, xpd = NA
)

