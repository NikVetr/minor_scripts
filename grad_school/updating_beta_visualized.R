#specify shape of beta prior
a_prior <- 1
b_prior <- 1

#count of successes and failures
n_success <- 0.65
n_fail <- 0.35

#compute shape of posterior distribution (also beta)
a_posterior <- a_prior + n_success
b_posterior <- b_prior + n_fail

#compute expected values
expected_value_prior <- a_prior / (a_prior + b_prior)
expected_value_posterior <- a_posterior / (a_posterior + b_posterior)

#turn it into a graph
prob_range <- 0:100 / 100
prior_dens <- dbeta(prob_range, a_prior, b_prior)
posterior_dens <- dbeta(prob_range, a_posterior, b_posterior)
plot(prob_range, prior_dens, type = "l", col = "purple", ylab = "probability density", xlab = "probability", 
     main = "Updating a Beta Distribution", ylim = c(0, max(c(prior_dens, posterior_dens))))
points(expected_value_prior, prior_dens[which.min(abs(prob_range - expected_value_prior))], col = "purple", pch = 4)
lines(prob_range, posterior_dens, col = "blue")
points(expected_value_posterior, posterior_dens[which.min(abs(prob_range - expected_value_posterior))], col = "blue", pch = 4)
abline(v = expected_value_prior, lty = 3, lwd = 0.5, col = "purple")
abline(v = expected_value_posterior, lty = 3, lwd = 0.5, col = "blue")
legend(x = "topright", legend = c("prior", "posterior", "expected value"), col = c("purple", "blue", "black"), lwd = c(1,1,0.5), lty = c(1,1,3))
text(labels = round(expected_value_prior, 2), expected_value_prior + 0.035, y = max(c(prior_dens, posterior_dens)), col = "purple")
text(labels = round(expected_value_posterior, 2), expected_value_posterior + 0.045, y = 0, col = "blue")

#visualize partial success uncertainty
f <- function(){
  success <- rbinom(1, 1, 0.65)
  rbeta(1, 1 + success, 2 - success)
}
mean(replicate(1E5, f()))
