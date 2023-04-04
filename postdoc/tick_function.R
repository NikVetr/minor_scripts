nice_ticks <- function(d, nt = 5){
  magvar <- 10^floor(log10(diff(range(d)) / nt))
  vals <- seq(from = min(d) - min(d) %% magvar, 
              to = max(d) - max(d) %% magvar, 
              by = magvar)
  tick_sub <- floor(length(vals) / ny)
  tick0 <-  which.min(abs(vals))
  return(vals[c(rev(seq(from = tick0, to = 1, by = -tick_sub)[-1]), tick0, seq(from = tick0, to = length(vals), by = tick_sub)[-1])])
}

nice_ticks(d = rnorm(100), 6)
