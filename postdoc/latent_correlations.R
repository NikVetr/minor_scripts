n <- 1E3
d <- 100

latent_cor <- function(n, d, es = 1){
  x <- rnorm(n)
  e <- matrix(rnorm(n*d), n, d) * es
  y <- x + e
  evls <- eigen(cov(y))
  v <- evls$vectors
  l <- evls$values
  cor(x, y %*% v[,1])
}



ds <- ceiling(1.5^(0:17))
cors <- abs(unlist(parallel::mclapply(ds, latent_cor, n = n, es = 2, mc.cores = 12)))
plot(ds, log10(1-cors), type = "l")
