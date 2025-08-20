#a script to test how to maximize diversity in sample selection

#functions
sample_block_R <- function(p, n_blocks, base_r, block_r) {
  
  # Set up random block sizes
  block_cuts <- abs(rnorm(n_blocks))^(1/2)
  block_cuts <- cumsum(ceiling(c(0, block_cuts / sum(block_cuts) * p)))
  block_cuts <- c(block_cuts[block_cuts < p], p)
  
  # Create base correlation matrix
  R <- diag(p) + base_r - diag(p) * base_r
  
  # Apply block-specific correlation structure
  for (i in 1:(length(block_cuts) - 1)) {
    block_size <- block_cuts[i + 1] - block_cuts[i]
    block_corr <- block_r[i]  # correlation for this block
    R[(block_cuts[i] + 1):block_cuts[i + 1], (block_cuts[i] + 1):block_cuts[i + 1]] <- 
      diag(block_size) + block_corr - diag(block_size) * block_corr
  }
  
  return(R)
}

utri <- function(x) x[upper.tri(x)]

#hyperparams
n <- 2E4
k <- length(LETTERS)
p <- 500
n_blocks <- 30
big_freq <- 0.95
subset_props <- c(big_freq, rep((1-big_freq) / (k-1), k-1))
subset_pops <- n * subset_props
group_labels <- rep(1:k, subset_pops)

#subset correlation matrices
within_Rs <- lapply(1:k, function(ki) sample_block_R(p = p, n_blocks = n_blocks, 
                                                     base_r = runif(1, 0, 0.2), 
                                                     block_r = runif(nblocks, 0.2, 0.9)))
within_Ls <- lapply(within_Rs, function(R) t(chol(R)))

#population means (iid standard normal)
pop_means <- t(replicate(k, rnorm(p)))

#sample data
x_list <- lapply(1:k, function(ki){
  pop_dev <- t(within_Ls[[ki]] %*% matrix(rnorm(p * subset_pops[ki]), ncol = subset_pops[ki], nrow = p))
  return(t(t(pop_dev) + pop_means[ki,]))
})
x <- do.call(rbind, x_list)
total_R <- cor(x)
# plot(utri(total_R), utri(within_Rs[[1]]))

#transform to PC space
PCA <- prcomp(x)
PCs <- PCA$x
PC_split <- split(PCs, group_labels)

#### plot PC scores ####
png(filename = "~/PC_stratification.png", 
                     width = 4096, height = 4096, pointsize = 25)
par(mfrow = c(10, 10), mar = c(4,4,2,2))

for(i in 1:100){
  cat(paste0(i, " "))
  plot(PCs[,i], group_labels,
       main = paste0("PC Axis ", i),
       ylab = "Group Index",
       xlab = "PC Score", yaxt = "n")
  axis(2, at = 1:k, labels = 1:k, las = 1)
  
}
dev.off()
