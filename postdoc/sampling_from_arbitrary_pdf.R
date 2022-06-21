f <- function(x) 9*x*exp(-3*x)
integrate(f, 0, Inf) #let's just check lol
q_f <- function(q) 1 / 3 * (-lamW::lambertWm1((q-1)/exp(1)) - 1) 

set.seed(1)
n <- 5E4
u <- runif(n, 0, 1)
x <- q_f(u)

hist(x, probability = T, breaks = 100)
curve(f, 0, 5, add = T, lwd = 2, col = "blue")