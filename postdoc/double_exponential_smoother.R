n <- 1E2
x <- 1:n
interpolate_at <- seq(min(x), max(x), length.out = 5E3)
y <- cumsum(rnorm(n, mean = 0, sd = sqrt(1E3 / n)))
dexp_smooth <- function(x, y, r, reweight_trunc_tail = F, reweight_n_pts = F, fix_endpoints = F, interpolate_at = NA){
  
  nx <- length(x)
  if(is.na(interpolate_at[1])){
    n <- nx
    w <- t(sapply(1:(n-1), function(i) c(dexp((x[i] - x[0:(i-1)]), rate = r), dexp(0, rate = r), dexp(-(x[i] - x[(i+1):n]), rate = r))))
    w <- rbind(w, dexp(x[n] - x, rate = r))
  } else {
    n <- length(interpolate_at)
    w <- t(sapply(1:(n-1), function(i) c(dexp((interpolate_at[i] - x[x<=interpolate_at[i]]), rate = r), 
                                         dexp((x[x>interpolate_at[i]] - interpolate_at[i]), rate = r))))
    w <- rbind(w, dexp(x[nx] - x, rate = r))
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(reweight_trunc_tail){
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:(n-1), function(i) c(pexp((x[i] - x[1]), rate = r), pexp(-(x[i] - x[n]), rate = r))))
    } else {
      tw <- t(sapply(1:(n-1), function(i) c(pexp((interpolate_at[i] - x[1]), rate = r), pexp(-(interpolate_at[i] - x[nx]), rate = r))))
    }
    tw[1,] <- tw[2,]
    tw <- rbind(tw, tw[n-1,])
    tw <- 1/tw
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], ri), rep(1, n-ri)) * c(rep(1, ri), rep(tw[ri,2], n-ri)))) #eh let's give it to the nth pt too
    } else {
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], sum(interpolate_at[ri] >= x)), rep(1, nx-sum(interpolate_at[ri] >= x))) * 
                                       c(rep(1, sum(interpolate_at[ri] >= x)), rep(tw[ri,2], nx-sum(interpolate_at[ri] >= x)))))
    }
    w <- t(sapply(1:n, function(ri) w[ri,] * tw[ri,]))
    
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(reweight_n_pts){
    if(is.na(interpolate_at[1])){
      tw <- cbind(0:(n-1), (n-1):0)
      tw <- 1/tw
    } else {
      tw <- sapply(1:n, function(i) sum(interpolate_at[i] >= x))
      tw <- cbind(tw, nx - tw)
    }
    mintw <- apply(tw, 1, min)
    tw[,1] <- tw[,1] / mintw 
    tw[,2] <- tw[,2] / mintw 
    
    tw[1,] <- tw[2,]; tw[n,] <- tw[n-1,]

    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], ri), rep(1, n-ri)) * c(rep(1, ri), rep(tw[ri,2], n-ri)))) #eh let's give it to the nth pt too
    } else {
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], sum(interpolate_at[ri] >= x)), rep(1, nx-sum(interpolate_at[ri] >= x))) * 
                       c(rep(1, sum(interpolate_at[ri] >= x)), rep(tw[ri,2], nx-sum(interpolate_at[ri] >= x)))))
    }
    
    w <- t(sapply(1:n, function(ri) w[ri,] * tw[ri,]))
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(fix_endpoints){
    w[1,] <- c(1, rep(0,nx-1))
    w[n,] <- c(rep(0,nx-1), 1)
  }
  
  wy <- c(w %*% t(t(y)))
  wy <- wy / apply(w,1,sum)
  return(wy)
}

par(mfrow = c(4,1))

r = 0.1
plot(x, y, type = "l")
# lines(x, dexp_smooth(x, y, 0.1), type = "l", col = 2, lwd = 2)
plot(x, y, type = "l", col = "grey75")
lines(x, dexp_smooth(x, y, r, F,F,F), type = "l", col = 2, lwd = 3)

plot(x, y, type = "l", col = "grey75")
lines(interpolate_at, dexp_smooth(x, y, r, F,F,F, interpolate_at = interpolate_at), type = "l", col = 2, lwd = 3)
plot(-diff(diff(dexp_smooth(x, y, r, interpolate_at = interpolate_at))), type = "l")
