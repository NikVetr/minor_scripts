n <- 20
loc <- orig_loc <- rnorm(n)
sep <- 0.025
pow <- 1
dists <- as.matrix(dist(loc))
overlap <- any(dists[upper.tri(dists)] < sep)
force <- 2
max_step <- 0.1

while(overlap){
  # print(loc)
  dirs <- (sapply(loc, function(loci) loci > loc) * 2 - 1)
  forces <- dists^(1/pow) * force
  forces[dists > (sep*2)] <- forces[dists > (sep*2)] / 10
  force_loc <-  apply(forces * dirs, 2, sum)
  maxf <- max(abs(force_loc))
  t <- sqrt((max_step * sep * 2) / maxf)
  disps <- force_loc * 0.5 * t^2
  loc <- loc + disps
  dists <- as.matrix(dist(loc))
  overlap <- any(dists[upper.tri(dists)] < sep)
}

plot(orig_loc, col = 3, pch = 19, x = rep(0, n), ylim = range(c(orig_loc, loc)), xlim = c(-1,2))
points(loc, col = 2, pch = 19, x = rep(1, n))
segments(0, y0 = orig_loc, 1, y1 = loc)
min(dists[upper.tri(dists)])
