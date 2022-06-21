ts <- c(1,2,4,8)
n <- 1E4
x1 <- t(apply(sapply(ts, function(t) rnorm(n, sd = sqrt(t))), 1, cumsum))
x2 <- t(apply(sapply(ts, function(t) rnorm(n, sd = sqrt(t))), 1, cumsum))
par(mfrow = c(3,1))
hist(sapply(1:n, function(i) cor(x1[i,], x2[i,])), main = "raw pearson corrs")
hist(sapply(1:n, function(i) cor(diff(c(0,x1[i,])) / sqrt(diff(c(0,ts))), diff(c(0,x2[i,])) / sqrt(diff(c(0,ts))))), 
     main = "differenced pearson corrs")

correlate_time_series <- function(x1, x2, relative_timepoints = c(0,1,2,4,8)){
  if(length(x1) != length(x2) || length(x1) != length(relative_timepoints)){
    stop("x1, x2, and relative_timepoints all have to be the same length")
  }
  time_diffs <- diff(relative_timepoints)
  x1_std_diffs <- diff(x1) / sqrt(time_diffs)
  x2_std_diffs <- diff(x2) / sqrt(time_diffs)
  return(cor(x1_std_diffs, x2_std_diffs))
}

hist(sapply(1:n, function(i) correlate_time_series(c(0,x1[i,]), c(0,x2[i,]), c(0,ts))),
     main = "differenced pearson corrs")
