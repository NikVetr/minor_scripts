#specify generative model parameters
n = 1E1
a1 = 3
b1 = 2
a2 = 5
b2 = 6

#simulate data
x1 = rnorm(n)
x2 = rnorm(n)
y1 = a1 + b1 * x1 + rnorm(n)
y2 = a2 + b2 * x2 + rnorm(n)
d <- data.frame(x = c(x1, x2), y = c(y1, y2), sex = c(rep(0, n), rep(1, n)))

#fit interaction model
m_inter <- lm(y ~ 1 + x + sex + x*sex, data = d)

#fit explicit separation model in single model
# lm(y ~ 0 + 1*sex + 1*(1-sex) + x*sex * x*(1-sex), data = d) #grr stop processing my formulae

#fit explicit separation model in single model, but avoiding lm()'s processing
d2 <- d
d2$xsex <- d2$x * d2$sex
d2$asex <- 1-d2$sex
d2$xasex <- d2$x * d2$asex
m_oneline <- lm(y ~ 0 + sex + asex + xsex + xasex, data = d2)

#fit "independent" / separate separation models 
m_sex0 <- lm(y ~ 1 + x, data = d[d$sex == 0,])
m_sex1 <- lm(y ~ 1 + x, data = d[d$sex == 1,])

#compare ols coefficient estimates
summary(m_inter)
summary(m_oneline)
summary(m_sex0)
summary(m_sex1)

eps = .Machine$double.eps ^ 0.5 #just to catch floating point precision

#compare concatenated and separate models
abs(m_oneline$coefficients["sex"] - m_sex1$coefficients["(Intercept)"]) < eps
abs(m_oneline$coefficients["asex"] - m_sex0$coefficients["(Intercept)"]) < eps
abs(m_oneline$coefficients["xsex"] - m_sex1$coefficients["x"]) < eps
abs(m_oneline$coefficients["xasex"] - m_sex0$coefficients["x"]) < eps

#compare concatenated and interaction model
abs(m_inter$coefficients["(Intercept)"] - m_oneline$coefficients["asex"]) < eps
abs(m_inter$coefficients["sex"] + m_inter$coefficients["(Intercept)"]  - m_oneline$coefficients["sex"]) < eps
abs(m_inter$coefficients["x"] - m_oneline$coefficients["xasex"]) < eps
abs(m_inter$coefficients["x"] + m_inter$coefficients["x:sex"]  - m_oneline$coefficients["xsex"]) < eps

#all true!

var(d2$xsex) 
var(d$x[d$sex == 1])

#### add a predictor to the dispersion term ####
library(dglm)

#specify generative model parameters
n = 1E3
a1 = 3
b1 = 2
a2 = 5
b2 = 6
s1 = 1
s2 = 2

#simulate data
x1 = rnorm(n)
x2 = rnorm(n)
y1 = a1 + b1 * x1 + rnorm(n, sd = s1)
y2 = a2 + b2 * x2 + rnorm(n, sd = s2)
d <- data.frame(x = c(x1, x2), y = c(y1, y2), sex = c(rep(0, n), rep(1, n)))

#fit model
m_dispersion <- dglm(y ~ x + sex + x*sex, dformula = ~ sex,  data = d)
summary(m_dispersion)
