
mu <- 0.5
vu <- 100
a <- 0 + mu * vu
b <- 0 + (1-mu) * vu
w <- 0.5
curve(expr = w * (dbeta(x, a, b)) + (1-w) * dbeta(x, b, a), 
      from = 0, to = 1, n = 2048, 
      main = paste0("a = ", a, ", b = ", b, ",\nmu = ", mu, ", vu = ", vu),
      ylab = "density")
curve(expr = w * (dbeta(x, a-30, b-30)) + (1-w) * dbeta(x, b-30, a-30), 
      from = 0, to = 1, n = 2048, 
      main = paste0("a = ", a, ", b = ", b, ",\nmu = ", mu, ", vu = ", vu),
      ylab = "density", col = 2, add = T)


# Define the components of the mixture
x <- seq(0, 1, length = 1000) # Values over which to evaluate the distributions

# Beta(15, 15)
beta1 <- dbeta(x, 15, 15)

# Beta(0.5, 0.5)
beta2 <- dbeta(x, 0.5, 0.5)
beta2[beta2 > 1E6] <- 2

# Uniform(0.499, 0.501)
uniform <- dunif(x, 0.495, 0.505)

# Assign equal weights to each component
w <- c(0.7, 0.25, 0.05)

# Mixture distribution
mixture <- w[1] * beta1 + w[2] * beta2 + w[3] * uniform

# Create the plot
plot(x, mixture, type = "n", col = "black", lwd = 2,
     ylab = "Density", xlab = "Probability", main = "Mixture Distribution")

# Draw the components using polygon()
polygon(c(x, rev(x)), c(w[1] * beta1, rep(0, length(beta1))),
        col = rgb(0, 0, 1, 0.4), border = NA) # ASE (Beta(15,15))
polygon(c(x, 1, rev(x), 0), c(w[2] * beta2, rep(0, length(beta2) + 1), w[2] * beta2[1]),
        col = rgb(1, 0.5, 0, 0.4), border = NA) # LoH (Beta(0.5,0.5))
polygon(c(x, rev(x)), c(w[3] * uniform, rep(0, length(uniform))),
        col = rgb(0, 1, 0, 0.4), border = NA) # AB (Uniform(0.499,0.501))

# Add the mixture curve
lines(x, mixture, col = "black", lwd = 2)

# Add a legend
legend("topright", legend = c("Mixture", "ASE", "LoH", "AB"),
       col = c("black", rgb(0, 0, 1, 0.4), rgb(1, 0.5, 0, 0.4), rgb(0, 1, 0, 0.4)),
       pch = 15, pt.cex = 2, bty = "n")

