params_1 <- c(1, 23)
nt <- 100

simulate_updating_betabinom <- function(params_1, nt){
  params <- matrix(0, nt, 2)
  obs <- rep(0, nt)
  params[1,] <- params_1
  for(i in 1:(nt-1)){
    theta[i] <- rbeta(1, params[i,1], params[i,2])
    obs[i] <- rbinom(1, 1, theta[i])
    params[i+1,] <- params[i,] + c(obs[i],1-obs[i])
  }
  return(list(theta = theta, obs = obs, params = params))
}

params_1[1] / sum(params_1)
mean(replicate(1E3, mean(simulate_updating_betabinom(params_1, nt)$obs)))

(1 - params_1[1] / sum(params_1))^nt
mean(replicate(1E3, all(simulate_updating_betabinom(params_1, nt)$obs == 0)))
