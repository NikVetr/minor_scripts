invlogit <- function(x){exp(x) / (1 + exp(x))}

mu <- 3
x <- rnorm(1E4, mu, sd = 1)
hist(x)

n <- 1E4
pres <- rbinom(n = n, size = 1, prob = invlogit((x-mu)*10+mu)) == 1
x_obs <- x[pres]
hist_obs <- hist(x_obs, breaks = 0:16/2, plot = F)
plot(hist_obs, xlim = c(0,8), ylim = c(0,2000), 
     main = paste0((1 - sum(pres) / length(pres))*100, "% of data is missing,\n", n - sum(pres), " missing observations / ", n, " total observations"), xlab = "observed values")
text(x = hist_obs$breaks + 0.25, y = hist_obs$counts + 50, labels = hist_obs$counts, cex = 0.75)
