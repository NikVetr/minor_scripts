n0 <- 1000

len0 <- 1/10/n0
pos <- runif(n0)
r <- 1E-4
t_n <- 1E4
coords <- cbind(pos - len0/2, pos + len0/2)

fuse <- function(coords){
  lens <- coords[,2] - coords[,1]
  cents <- (coords[,1] + coords[,2])/2
  ord_cents <- order(cents)
  cent_dists <- diff(cents[ord_cents])
  comb_lens <- (lens[ord_cents][1:(nrow(coords)-1)] + lens[ord_cents][2:(nrow(coords))]) / 2 
  to_fuse <- cent_dists < comb_lens
  rle2f <- rle(to_fuse)
  inds <- cumsum(c(0, rle2f$lengths))
  inds[1] <- 1
  fuse_inds <- which(rle2f$values)
  to_fuse_inds <- cbind(inds[fuse_inds], inds[fuse_inds+1]) + 1
}

for(i in 1:t_n){
  prop_disp 
}
