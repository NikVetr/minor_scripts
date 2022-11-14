mu <- 1.15^(-5:75)

n <- 10000
base_size <- 50
kappa <- 1
gamma <- 0
d <- sapply(mu, function(mu_i) rnbinom(n = n, mu = mu_i, size = kappa * mu_i^gamma * base_size))
# d <- sapply(mu, function(mu_i) rnbinom(n = n, mu = mu_i, size = base_size))
# 
# plot((t(apply(d, 2, quantile, c(0.01, 0.99))) / apply(d, 2, mean)))


par(mfrow = c(3,1))
plot(log(mu), log(apply(d, 2, var)), 
     main = latex2exp::TeX(paste0("$\\kappa = ", kappa, ", \\gamma = ", gamma, ", base\\_size = ", base_size, "$"))); abline(0,1)
plot(rep(mu, each = n), d)
plot(rep(log(mu), each = n), log(d))
