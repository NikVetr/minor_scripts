#meta-params
set.seed(1)
n <- 200
p <- 1000
np <- n * p

#simulate rin-independent variation
mus <- exp(rnorm(p, 1, 0.5))
sds <- exp(rnorm(p, 0.5, 0.25))
true_vals <- sapply(1:p, function(i) (rnorm(n, mus[i], sds[i])))

#simulate rin values
batch_counts <- ceiling(gtools::rdirichlet(1, rep(1, 20)) * np)
batch_counts <- batch_counts[cumsum(batch_counts) < np]
batch_counts <- c(batch_counts, np - sum(batch_counts))
rin <- exp(matrix(rnorm(np, 2, 0.1), n, p) + 
             matrix(rnorm(p, sd = 0.05), n, p, byrow = T) +
             matrix(rnorm(n, sd = 0.05), n, p, byrow = F) +
             matrix(rep(rnorm(length(batch_counts), sd = 0.05), batch_counts, each = T), n, p, byrow = T))

#simulate rin-dependent variation and add to rin-independent variation
rin_err <- matrix(rnorm(np, sd = 10 * (rin - min(rin))), n, p, byrow = F)
obs_vals <- (true_vals) + (rin - min(rin)) * 5 + rin_err
obs_vals <- apply(obs_vals, 2, scale)

#plot first two PCs
pca <- prcomp(obs_vals)
nlev <- 5
cols <- viridisLite::viridis(nlev)
mrin <- apply(rin, 1, mean)
rinval <- floor(rank(mrin) / length(mrin) * (nlev) - 1E-6) + 1
par(mfrow = c(2,2))
plot(pca$x[,1:2], col = adjustcolor(cols[rinval], 0.5), pch = 19)
legend(pch = 19, col = cols, lty = 1, x = "topright", 
       legend = sapply(1:nlev, function(i) paste0("(", 
                paste0(round(range(mrin[rinval == i]), 2), collapse = ", "), ")")),
       title = "mean rin value", cex = 0.5)

#plot bounding contours
kdes <- lapply(sort(unique(rinval)), function(i) 
  MASS::kde2d(pca$x[rinval == i,1], pca$x[rinval == i, 2]))
for(i in 1:nlev){contour(kdes[[i]], add = T, lwd = 2, col = cols[i], nlevels = 3, )}

#plot rin vs variance in sliding window
window_width <- diff(quantile(rin, c(0.4, 0.6)))
nwind <- 100
windows <- seq(from = min(rin), to = max(rin) - window_width, length.out = nwind)
windows <- data.frame(cbind(L = windows, R = windows + window_width))

#slide window to get sample expectations and variances
sample_vals <- parallel::mclapply(1:nwind, function(i){
  hits <- data.frame(which(rin > windows$L[i] & rin < windows$R[i], arr.ind = T))
  thits <- table(hits$col) 
  subhits <- hits[hits$col %in% as.integer(names(thits[thits > 5])),]
  subhits <- split(subhits, subhits$col)
  if(length(subhits) == 0){return(numeric(0))}
  cbind(t(sapply(subhits, function(x){
    vals <- apply(x, 1, function(j) obs_vals[j[1], j[2]])
    c(expectation = mean(vals), variance = var(vals) * length(vals) / (length(vals) - 1))
    })), bin = i)
}, mc.cores = 12)
sample_vals <- data.frame(do.call(rbind, sample_vals))

#plot these
plot(apply(windows, 1, mean)[sample_vals$bin], sample_vals$expectation, col = adjustcolor(1, 0.01), pch = 19,
     xlab = "rin window center", ylab = "expectation")
plot(apply(windows, 1, mean)[sample_vals$bin], sample_vals$variance, col = adjustcolor(1, 0.01), pch = 19,
     xlab = "rin window center", ylab = "variance")
