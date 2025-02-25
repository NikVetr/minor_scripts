library(copula)

ev_copula <- evCopula(family = "tev", param = 0.7)  # Example: Gumbel copula
n <- 100
u <- rCopula(n, ev_copula)  # Simulate uniform data

# Transform to desired marginal distributions (e.g., Exponential)
lambda <- 1000
data <- qpois(u, lambda = lambda)

#plot both
par(mfrow=c(1,2))
plot(u)
plot(log(data))

# fit extreme value copula

# Fit Gumbel copula as an extreme value copula
fit <- fitCopula(evCopula(family = "tev"), pobs(data), method = "ml")

# Apply the evTestA to test if data exhibits extreme value dependence
evTestA(pobs(data), N = 1000, derivatives = "Cn")
cor.test(data[,1], data[,2])[c("estimate", "p.value")]
