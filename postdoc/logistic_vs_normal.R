#functions
logit <- function(p) log(p/(1-p))
invlogit <- function(x) {exp(x) / (1+exp(x))}

#simulate normal outcomes w/ a binary predictor
n = 1E2
probs <- rbeta(n,2,2)
trait <- sapply(1:n, function(i) sample(c(0,1), size = 1, prob = c(probs[i], 1-probs[i])))
a = 4
b = 3
sigma = 1
age <- a + trait*b + rnorm(n, sd = sigma)

#symmetry here
summary(lm(age~trait))$coefficients
summary(lm(trait~age))$coefficients

#but not here
summary(glm(trait~age, family = "binomial"))$coefficients
summary(glm(trait~age, family = binomial(link = "probit")))$coefficients

#how's it looks if we model the logit probs?
logit_probs <- logit(probs)
summary(lm(logit_probs~age))$coefficients


#let's go the other way around and simulate bernoulli outcomes with a continuous predictor
a = 1
b = 3
x <- runif(n, -1, 1)
lp <- a + b*x + rnorm(n)
p <- invlogit(lp)
bern <- sapply(1:n, function(i) sample(c(0,1), size = 1, prob = c(p[i], 1-p[i])))

summary(lm(x~bern))$coefficients
summary(glm(bern~x, family = "binomial"))$coefficients

# summary(lm(bern~x))$coefficients
summary(lm(lp~x))$coefficients
