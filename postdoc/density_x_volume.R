# ks <- round((2^(1/2))^(1:40))
# stopping_dist <- sapply(ks, function(k) {
#   x <- -400:400/100
#   d <- dnorm(x)
#   dxv <- data.frame(t(sapply(0:((length(x)-1)/2-1), function(i){
#     ivals <- c((length(x)-1)/2-i, (length(x)-1)/2+i);
#     vol <- diff(x[ivals])^k
#     min_dens <- d[ivals[1]]^k
#     mass_contained <- (mean(d[ivals[1]:ivals[2]]) * diff(x[ivals]))^k
#     dxv <- vol * min_dens
#     dists <- sqrt(x[ivals[1]]^2*k)
#     c(vol = vol, min_dens = min_dens, dxv = dxv, dists = dists, mass = mass_contained)
#   })))
#   dxv$dists[which.max(dxv$dxv)]
# })
# plot(ks, stopping_dist)


#now on log scale
x <- -500:500/100
ks <- round((2^(1/2))^(1:40))

#E(x) = 0, var(x) = 1
ld <- dnorm(x, log = T)
ld <- dexp(abs(x), log = T) #laplace / gumbel
ld <- dunif(x, min = -sqrt(3), max = sqrt(3), log = T)
ld <- dlogis(x, scale = sqrt(3/pi^2), log = T)

#E(x) = 0
ld <- dt(x, df = 3, log = T)
ld <- dcauchy(x, location = 0, scale = 1, log = T)

stopping_dist <- sapply(ks, function(k) {
  ldxv <- data.frame(t(sapply(0:((length(x)-1)/2-1), function(i){
    ivals <- c((length(x)-1)/2-i, (length(x)-1)/2+i);
    lvol <- log(diff(x[ivals]))*k
    min_ldens <- ld[ivals[1]]*k
    lmass_contained <- (log(mean(exp(ld[ivals[1]:ivals[2]]))) + diff(x[ivals])) * k
    ldxv <- lvol + min_ldens
    ldists <- (log(abs(x[ivals[1]])) * 2 + log(k)) / 2
    c(lvol = lvol, min_ldens = min_ldens, ldxv = ldxv, ldists = ldists, lmass = lmass_contained)
  })))
  ldxv$ldists[which.max(ldxv$ldxv)]
})

par(mar = c(6,6,3,2))
plot(log10(ks), log10(exp(stopping_dist)), type = "l", xlab = latex2exp::TeX("k (dimension)"), xaxt = "n", 
     ylab = "optimal stopping distance (maximizing minimum density x volume)",
     yaxt = "n", lwd = 2, 
     main = paste0(c("intercept = ", "slope = "), 
                   round(lm(log10(exp(stopping_dist)) ~ log10(ks))$coefficients, 2), 
                   collapse = ", ")
     )
xt <- 0:floor(log10(max(ks)))
axis(1, xt, labels = latex2exp::TeX(paste0("10$^{", xt, "}$")))
yt <- seq(0, max(log10(exp(stopping_dist))), by = 0.5)
axis(2, yt, labels = latex2exp::TeX(paste0("10$^{", yt, "}$")))
abline(h = yt, lty = 3, cex = 0.5, col = adjustcolor(1,0.5))
abline(v = xt, lty = 3, cex = 0.5, col = adjustcolor(1,0.5))
