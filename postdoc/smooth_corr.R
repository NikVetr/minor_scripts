R <- sample.R

smoooth_EVs <- function(R, n){
  evs <- eigen(R)
  p <- dim(R)[1]
  logEVs <- log(evs$values[1:(n-1)])
  spline_fit <- smooth.spline(1:(n-1), logEVs)
  logEVs <- c(logEVs, predict(spline_fit, n:p)$y)
  newEVs <- exp(logEVs)
  newEVs <- newEVs / sum(newEVs) * p
  newR <- (evs$vectors) %*% diag(newEVs) %*% t(evs$vectors)
  return(cov2cor(newR))
}
