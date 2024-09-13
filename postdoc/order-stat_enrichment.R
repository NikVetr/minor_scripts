# set.seed(1)
p <- 30 #analogous to cell types
n <- 5000 #analogous to genes or analytes

#all genes x cell types have the same true effect (0)
#but these effects are observed with noise (aka estimated)
d <- matrix(rnorm(p*n), nrow = n, ncol = p)

#correlated ps?
base_r <- 0.1
block_r <- 0.7
n_blocks <- 10 #approximately
#random block sizes (roughly p * dirichlet(2,2,2...) distributed)
block_cuts <- abs(rnorm(n_blocks))^(1/2)
block_cuts <- cumsum(ceiling(c(0,block_cuts / sum(block_cuts) * p)))
block_cuts <- c(block_cuts[block_cuts < p], p)
R <- diag(p) + base_r - diag(p) * base_r
for(i in 1:(length(block_cuts)-1)){
  block_size <- block_cuts[i+1] - block_cuts[i]
  R[(block_cuts[i]+1):block_cuts[i+1], (block_cuts[i]+1):block_cuts[i+1]] <- diag(block_size) + block_r - diag(block_size) * block_r
}
d <- d %*% chol(R)

#we threshold on the magnitude of the effect
thresh <- 2
d_nthresh <- apply(d, 1, function(x) sum(abs(x) > thresh))

#and find the max effects in each gene, splitting by post-threshold set size
d_max <- apply(d, 1, function(x) max(abs(x)))
d_mean_pass <- apply(d, 1, function(x) mean(abs(x[x>thresh])))
d_mean_pass_thresh <- split(d_mean_pass[d_nthresh>0], d_nthresh[d_nthresh>0])
d_mean_pass_thresh <- d_max_thresh[sapply(d_mean_pass_thresh, length) > 1] #remove sets with only one gene
d_max_thresh <- split(d_max, d_nthresh)
d_max_thresh <- d_max_thresh[sapply(d_max_thresh, length) > 1] #remove sets with only one gene

#what's the mean effect across max effects? and the SE of the mean + confidence interval?
mean_d_max_thresh <- sapply(d_max_thresh, mean)
SE_mean_d_max_thresh <- sapply(d_max_thresh, function(x) sd(x)/sqrt(length(x)))
CI_mean_d_max_thresh <- rbind(lower = mean_d_max_thresh - 2 * SE_mean_d_max_thresh,
                              upper = mean_d_max_thresh + 2 * SE_mean_d_max_thresh)

mean_d_mean_thresh <- sapply(d_mean_pass_thresh, mean)
SE_mean_d_mean_thresh <- sapply(d_mean_pass_thresh, function(x) sd(x)/sqrt(length(x)))
CI_mean_d_mean_thresh <- rbind(lower = mean_d_mean_thresh - 2 * SE_mean_d_mean_thresh,
                              upper = mean_d_mean_thresh + 2 * SE_mean_d_mean_thresh)

#do some plotting
par(mar = c(5,7,2,2), mfrow = c(2,1))
nthresh <- as.integer(names(mean_d_max_thresh))
plot(nthresh, mean_d_max_thresh, ylim = range(CI_mean_d_max_thresh),
     xlab = "number of samples above threshold", 
     ylab = "mean of max across entire sample\n(both above and below threshold)")
segments(x0 = nthresh, y0 = CI_mean_d_max_thresh[1,], y1 = CI_mean_d_max_thresh[2,])
legend("bottomright", legend = c("mean", "mean +/- 2SE"), lwd = c(NA,1), pch = c(1,NA))

nthresh <- as.integer(names(mean_d_mean_thresh))
plot(nthresh, mean_d_mean_thresh, ylim = range(CI_mean_d_mean_thresh),
     xlab = "number of samples above threshold", 
     ylab = "mean of mean above threshold")
segments(x0 = nthresh, y0 = CI_mean_d_mean_thresh[1,], y1 = CI_mean_d_mean_thresh[2,])
legend("bottomright", legend = c("mean", "mean +/- 2SE"), lwd = c(NA,1), pch = c(1,NA))
