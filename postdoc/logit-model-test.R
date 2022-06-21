logit <- function(p) log(p/(1-p))
invlogit <- function(x) exp(x) / (1+exp(x))

n <- 2E3
sex <- rbinom(n, 1, 0.5)
stem <- rbinom(n, 1, invlogit(-1.5 + sex*3))

stem_b <- 1.5
sex_b <- 0.75
stem_x_sex_b <- 1
intercept <- -1.25

success <- rbinom(n, 1, invlogit(intercept + sex_b * sex + stem_b * stem + stem * sex * stem_x_sex_b))

glm(formula = success ~ 1 + stem + sex + stem * sex, family = binomial(link = "logit"))
