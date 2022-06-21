sim_marginal_ll <- function(){

  prior_sd = 2
  data_sd = 3
  
  x <- rnorm(n = 1, mean = 1, sd = 4)
  marginal_x_dens_analytic <- dnorm(x, mean = 0, sd = sqrt(prior_sd^2 + data_sd^2), log = T) #compute marginal density of x
  
  
  xrange <- -2000:2000/ 100
  marginal_x_dens_numeric <- log(sum(sapply(xrange, function(xi) dnorm(x, mean = xi, sd = data_sd, log = F) #compute density of x in data dist
                                                                 * dnorm(xi, mean = 0, sd = prior_sd, log = F))) # weigh by  density of mean of data dist in prior dist
                                                                 / sum(dnorm(xrange, 0, prior_sd, F))) #rescale so things sum to 1
  
  c(numeric = marginal_x_dens_numeric, analytic = marginal_x_dens_analytic)


}

#behold! identitical
plot(t(replicate(20, sim_marginal_ll())), cex = 2, pch = 19, col = adjustcolor(1, 0.5), main = "marginal likelihood")
abline(0,1, col = 2, lwd = 3)
