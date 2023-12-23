#libraries
library(transport)
library(philentropy)
library(optimx)

#### functions ####
trim <- function(x, n=1) x[(n+1):(length(x)-n)]
xentropy <- function(x, y, logy = F){
  ifelse(logy, -sum(x*y), -sum(x*log(y)))
}
moments_from_grid <- function(x, d, m = c(1:2), ll = F){
  
  #prevent underflow
  if(ll){
    d <- d - max(d) + 10
    d <- exp(d)
  } 
  
  #compute first and later moments
  m1 <- sum(x*d) / sum(d)
  m_eval <- setdiff(m, 1)
  ms <- sapply(m_eval, function(m_i) sum((x-m1)^m_i*d) / sum(d))
  if(1 %in% m){return(unlist(c(m1, ms)))}else{return(ms)}  
  
}

#test function internals betw two dists
max_desired_mass <- 0.01
quantiles <- seq(max_desired_mass,1-max_desired_mass,by=max_desired_mass)
xq <- qlnorm(quantiles)
yq <- qgamma(quantiles, shape = 1)

# xq <- qnorm(quantiles)
# yq <- qnorm(quantiles)



rxy <- range(c(xq, yq))
min_interval <- min(c(diff(xq), diff(yq)))
jvals <- c(-Inf, seq(rxy[1], rxy[2], by = min_interval), Inf)

xm <- diff(plnorm(jvals))
ym <- diff(pgamma(jvals, shape = 1))

# xm <- diff(pnorm(jvals))
# ym <- diff(pnorm(jvals))


max_obs_mass <- max(c(xm, ym))
while(max_obs_mass > max_desired_mass){
  min_interval <- min_interval * 0.99
  jvals <- c(-Inf, seq(rxy[1], rxy[2], by = min_interval), Inf)
  xm <- diff(plnorm(jvals))
  ym <- diff(pgamma(jvals, shape = 1))
  
  # xm <- diff(pnorm(jvals))
  # ym <- diff(pnorm(jvals))
  
  max_obs_mass <- max(c(trim(xm), trim(ym)))
}



c(
  wasserstein = wasserstein1d(a = xm, b = ym),
  KL = KL(x = rbind(xm, ym), test.na = F, unit = "log2", epsilon = 1E-6),
  JS = JSD(x = rbind(xm, ym), test.na = F, unit = "log2"),
  bhattacharyya = bhattacharyya(P = xm, Q = ym, testNA = F, unit = "log2", epsilon = 1E-6),
  xentropy = xentropy(x = xm, y = ym)
)

#can also use method of moments? find the moments of one distribution, set the other to it, solve for parameters
#or maximum likelihood? find quantiles of starting dist, weigh them by their density in target dist, optimize params of target
#hm this might be equiv to using x-entropy actually


#find optimal gamma for log(std-normal)
prob_dist <- function(par, data){

  #extract parameters
  metric <- data[["metric"]]
  max_desired_mass <- data[["max_desired_mass"]]
  
  dist_1 <- data[["dist_1"]]
  dist_2 <- data[["dist_2"]]
  dist_1_pars <- data[["dist_1_pars"]]
  dist_2_pars <- c(unlist(data[["dist_2_pars"]]), unlist(par))
  
  #evaluate quantiles
  quantiles <- seq(max_desired_mass,1-max_desired_mass,by=max_desired_mass)
  xq <- do.call(paste0("q", dist_1), c(list(p = quantiles), as.list(dist_1_pars)))
  yq <- do.call(paste0("q", dist_2), c(list(p = quantiles), as.list(dist_2_pars)))
  rxy <- range(c(xq, yq))
  min_interval <- min(c(diff(xq), diff(yq)))
  jvals <- c(-Inf, seq(rxy[1], rxy[2], by = min_interval), Inf)
  
  #evaluate probability masses
  xm <- diff(do.call(paste0("p", dist_1), c(list(q = jvals), as.list(dist_1_pars))))
  ym <- diff(do.call(paste0("p", dist_2), c(list(q = jvals), as.list(dist_2_pars))))
  
  #get quantile -> bins below the desired threshoild
  
  max_obs_mass <- max(c(xm, ym))
  while(max_obs_mass > max_desired_mass){
    min_interval <- min_interval * 0.99
    jvals <- c(-Inf, seq(rxy[1], rxy[2], by = min_interval), Inf)
    xm <- diff(do.call(paste0("p", dist_1), c(list(q = jvals), as.list(dist_1_pars))))
    ym <- diff(do.call(paste0("p", dist_2), c(list(q = jvals), as.list(dist_2_pars))))
    max_obs_mass <- max(c(trim(xm), trim(ym)))
  }
  
  #evaluate divergence / distance metric
  if(metric == "KL"){
    out <- suppressMessages(KL(x = rbind(xm, ym), test.na = F, unit = "log2", epsilon = 1E-6))
  } else if(metric == "wasserstein"){
    out <- wasserstein1d(a = xm, b = ym) * min_interval * 1E5
  } else if(metric == "JS"){
    out <- suppressMessages(JSD(x = rbind(xm, ym) , test.na = F, unit = "log2"))
  } else if(metric == "bhattacharyya"){
    out <- bhattacharyya(P = xm, Q = ym, testNA = F, unit = "log2", epsilon = 1E-6)
  } else if(metric == "xentropy"){
    out <- xentropy(x = xm, y = ym + 1E-6, logy = F)
  }
  
  print(paste0(c(paste0(names(par), " = ", par), paste0(metric, " = ", round(out, 10))), collapse = ", "))
  return(out * length(jvals))
  
}

data <- list(metric = "wasserstein",
             max_desired_mass = 0.01,
             dist_1 = "lnorm",
             dist_2 = "gamma",
             dist_1_pars = c(meanlog = 5, sdlog = 0.5),
             dist_2_pars = NULL)
par <- c(shape = 2, rate = 2)

#find moments
#lnorm
# mean(rlnorm(1E5, dist_1_pars["meanlog"], dist_1_pars["sdlog"]))
lnorm_mean <- exp(data$dist_1_pars["meanlog"] + data$dist_1_pars["sdlog"]^2 / 2)

# var(rlnorm(1E6, dist_1_pars["meanlog"], dist_1_pars["sdlog"]))
lnorm_var <- exp(2 * data$dist_1_pars["meanlog"] + data$dist_1_pars["sdlog"]^2) * (exp(data$dist_1_pars["sdlog"]^2) - 1)

#gamma
par <- c(shape = as.numeric(lnorm_mean^2 / lnorm_var), rate = as.numeric(lnorm_mean / lnorm_var))
lnorm_mean
par["shape"] / par["rate"] #mean
lnorm_var
par["shape"] / par["rate"]^2 #variance


fit <- optimx::optimx(par = par,
      fn = prob_dist,
      data = data,
      method = "L-BFGS-B", 
      upper = c(Inf, Inf), 
      lower = c(0, 0)
)

xr <- 1:10000/10
plot(xr, dlnorm(xr, data$dist_1_pars["meanlog"], data$dist_1_pars["sdlog"]), type = "l")
lines(xr, dgamma(xr, shape = par["shape"], rate = par["rate"]), col = 3)
lines(xr, dgamma(xr, shape = fit$shape, rate = fit$rate), col = 2)
lines(xr, dgamma(xr, shape = 1.005, rate = 3.86), col = 2)
